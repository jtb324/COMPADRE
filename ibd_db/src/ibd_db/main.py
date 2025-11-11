from fastapi import FastAPI, Depends, HTTPException, Request, status
import redis.asyncio as aredis
import redis
from redis.backoff import ExponentialBackoff
from redis.retry import Retry
from contextlib import asynccontextmanager


retry = Retry(ExponentialBackoff(), 8)


# The async context manager function as our connection pool that is made once for the lifespan of the application
@asynccontextmanager
async def lifespan(app: FastAPI):
    print("initializing redis pool")
    app.state.redis = await aredis.from_url(
        "redis://localhost:6379", decode_responses=True, retry=retry
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
