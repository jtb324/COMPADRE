from fastapi import FastAPI
from fastapi.responses import JSONResponse

app = FastAPI()


@app.get("/")
def read_root():
    return {"Success": "Hello World"}


@app.get("/ersa/{data}")
def run_ersa():
    return {"Success": "Ran ERSA"}


@app.get("/pca/")
def run_pop_classifier():
    return {"Success": "Ran population classifier"}


@app.get("/padre/")
def run_padre():
    return {"Success": "Ran padre"}
