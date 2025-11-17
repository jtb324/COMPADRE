from fastapi import FastAPI, Depends, HTTPException, Request, status
from fastapi.responses import StreamingResponse
import json
import redis.asyncio as aredis
import redis
from redis.backoff import ExponentialBackoff
from redis.retry import Retry
from contextlib import asynccontextmanager
import os

REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379")

retry = Retry(ExponentialBackoff(), 8)


# The async context manager function as our connection pool that is made once for the lifespan of the application
@asynccontextmanager
async def lifespan(app: FastAPI):
    print("initializing redis pool")
    app.state.redis = await aredis.from_url(
        REDIS_URL, decode_responses=True, retry=retry
    )
    try:
        yield
    finally:
        print("shutting down server and closing the redis pool")
        await app.state.redis.aclose()


# We can return a connection from the pool rather than remake it
async def get_redis_connection(request: Request) -> aredis.Redis:
    return request.app.state.redis


app = FastAPI(lifespan=lifespan)


@app.get("/")
async def read_root():
    return {"return_code": 200, "message": "FastAPI IBD segment server"}


@app.get("/health/")
async def health_check(redis_conn: aredis.Redis = Depends(get_redis_connection)):
    try:
        await redis_conn.ping()

        return {"status_code": status.HTTP_200_OK, "message": "Redis service is active"}
    except redis.exceptions.ConnectionError as e:
        return {
            "status_code": status.HTTP_503_SERVICE_UNAVAILABLE,
            "message": "Unable to reach the redis server",
        }


@app.get("/pair_count/")
async def pair_count(redis_conn: aredis.Redis = Depends(get_redis_connection)):

    try:
        dbsize = await redis_conn.dbsize()

        return {"status_code": status.HTTP_200_OK, "message": dbsize}
    except redis.exceptions.ConnectionError:
        return {
            "status_code": status.HTTP_503_SERVICE_UNAVAILABLE,
            "message": "Unable to reach the redis server",
        }


async def key_scanner(request: Request, batch_size: int):
    """
    Async generator function to scan keys and yield batches.
    """
    batch = []
    try:
        # Use scan_iter to create an async iterator
        async for key in request.app.state.redis.scan_iter(count=batch_size):
            value = await request.app.state.redis.lrange(key, 0, -1)
            batch.append({"pair_id": key, "segment_info": value})
            if len(batch) >= batch_size:
                # Yield the batch as a JSON line
                yield f"{json.dumps(batch)}\n"
                batch = []

        # Yield any remaining keys in the last batch
        if batch:
            yield f"{json.dumps(batch)}\n"

    except Exception as e:
        print(f"Error during scan: {e}")
        yield f'{{"error": "An error occurred during the scan: {e}"}}\n'


@app.get("/iterate_pairs")
async def scan_keys(request: Request, batch_size: int = 100):
    """
    Scans for keys in Redis and streams the results in batches.

    - `match`: The pattern to match (e.g., "user:*")
    - `batch_size`: The number of keys to return in each JSON array.
    """
    # Use StreamingResponse to send the data as it's generated
    return StreamingResponse(
        key_scanner(request, batch_size), media_type="application/x-json-stream"
    )
