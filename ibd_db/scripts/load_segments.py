import argparse
import sys
from xopen import xopen
from pathlib import Path
import redis
import json
from redis.retry import Retry
from redis.backoff import ExponentialBackoff
from typing import Any, Generator
from collections import namedtuple

SegmentInfo = namedtuple("SegmentInfo", ["id1", "id2", "values"])


# This object is going to help us keep track of whether or
# not the file has an ibd2 statue and what the minimum segment
# threshold is
class RuntimeState:
    def __init__(self, min_cm_threshold: float) -> None:
        self.ibd2_status = "false"
        self.min_cm_threshold = min_cm_threshold
        self.min_cm_encountered = -1.0


def convert_value(value):
    """Convert string values of ERSA options to appropriate Python types"""
    # Try int first
    try:
        return int(value)
    except ValueError:
        pass

    # Try float
    try:
        return float(value)
    except ValueError:
        pass

    # Handle empty string as None for file paths
    if value == "None":
        return None

    # Return as string
    return value


def clean_ersa_options(ersa_flag_options: str) -> dict[str, Any]:
    """check if any of the ersa_options are invalid and correct
    them if they are

    Parameters
    ----------
    ersa_flag_options : str
        string that has all of the ersa options for the analysis.
        Each command option is separated by a pipe '|'

    Returns
    -------
    dict[str, Any]
        returns a dictionary with sanitized correct values for ERSA
    """

    # Check data types and revert to default values if they're not correct
    additional_options = {}

    # First step -- read in ERSA options from COMPADRE runtime flags
    ersa_flags = ersa_flag_options.strip('"').split("|")
    flag_names = ersa_flags[::2]
    flag_values = ersa_flags[1::2]

    for flag, value in zip(flag_names, flag_values):
        additional_options[flag] = convert_value(value)

    # additional_options = clean_ersa_options(additional_options)

    int_options = {
        "ascertained_position": -1,
        "mask_region_sim_count": 0,
        "mask_region_cross_length": 1000000,
        "mask_region_threshold": 4,
        "number_of_chromosomes": 22,
        "max_meioses": 40,
        "parent_offspring_zscore": 2.33,
        "pois_mean": 13.73,
        "exp_mean": 3.197036753,
        "rec_per_meioses": 35.2548101,
        "min_cm": 2.5,
        "confidence_level": 0.9,
        "max_cm": 10,
    }

    sanitized_dict = {}

    for option, value in additional_options.items():
        if isinstance(value, (int, float)):
            sanitized_dict[option] = value
        else:
            new_val = int_options.get(option, value)
            print(
                f"WARNING: the option:value combo {option}:{value} is not valid. The value should be either a integer or a float. Replacing old value with the new value, {new_val}"
            )
            sanitized_dict[option] = new_val

    return sanitized_dict


def read_segment_data(
    segment_data_file: Path, runtime_state: RuntimeState
) -> Generator[SegmentInfo, None, None]:
    """function to read through the shared segment file and format
    each line.

    Parameter
    ---------
    segment_data_file : Path
        filepath to a tab separated text file that contains the pairwise IBD segments. This file can either be compressed with gzip or decompressed

    ersa_flag_str : str
        insert description

    Returns
    -------
    Iterator[str]
        yields each line in the file"""
    with xopen(
        segment_data_file, "r"
    ) as f:  # Open the large file here and populate dictionary that stays in system memory

        # Check number of columns to determine if it has IBD1/2 data or not

        for line in f:
            if "chrom" in line:  # skip header if it exists
                continue
            ls = line.strip().split("\t")

            # If the segment is shorter than the minimum centimorgan
            # threshold then we can just skip the line
            if float(ls[4]) <= runtime_state.min_cm_threshold:
                continue

            if len(ls) == 7:  # File with segment-by-segment IBD status in 7th column
                runtime_state.ibd2_status = "true"

                iid1, iid2, start, end, cmlen, chrom, ibd = (
                    ls[0],
                    ls[1],
                    int(ls[2]),
                    int(ls[3]),
                    round(float(ls[4]), 2),
                    int(ls[5]),
                    int(ls[6].strip()),
                )

            # This block is hit if we are dealing with GERMLINE2 data with
            # out using the haploid flag
            elif len(ls) == 6:  # File without segment-by-segment IBD status
                iid1, iid2, start, end, cmlen, chrom = (
                    ls[0],
                    ls[1],
                    int(ls[2]),
                    int(ls[3]),
                    round(float(ls[4]), 2),
                    int(ls[5].strip()),
                )
                ibd = "NA"

            # This block is hit if we are dealing with GERMLINE 1 data
            elif len(ls) == 11:
                iid1, iid2, start, end, cmlen, chrom = (
                    ls[0].split(" ")[0],
                    ls[1].split(" ")[0],
                    int(ls[3].split(" ")[0]),
                    int(ls[3].split(" ")[1]),
                    round(float(ls[6]), 2),
                    int(ls[2]),
                )

                ibd = "NA"
            else:
                print(
                    "ERROR: [COMPADRE] Unrecognized segment file format. Please refer to the README (https://github.com/belowlab/compadre) for formatting guidelines."
                )
                sys.exit(1)
            # lets store all of this values in a tuple
            value = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "cM": cmlen,
                "ibd": ibd,
            }  # no ibd1/2 data in this input

            # we need to record what the lowest cm encounter is because we will use that to adjust the state later
            if runtime_state.min_cm_encountered == -1:
                runtime_state.min_cm_encountered = cmlen
            elif runtime_state.min_cm_encountered > cmlen:
                runtime_state.min_cm_encountered = cmlen
            yield SegmentInfo(*[iid1, iid2, value])  # TODO: fix the type annotation


# Figure out the true min_cm to be passed into ERSA based on the contents of segment_dict
def add_values_to_db(info: SegmentInfo, db_handler: redis.Redis) -> None:
    """add the value of interest to the database"""
    pair_id = f"pair:{info.id1}:{info.id2}"
    rev_pair_id = f"pair:{info.id2}:{info.id1}"

    if db_handler.exists(pair_id):
        lookup_id = pair_id
    elif db_handler.exists(rev_pair_id):
        lookup_id = rev_pair_id
    else:
        lookup_id = pair_id

    serialized_data = json.dumps(info.values)

    try:
        db_handler.rpush(lookup_id, serialized_data)
    except redis.exceptions.ResponseError:
        print(
            f"Failed to add the value to the database. Error occurred for pair_id: {pair_id} - value {info.values}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        "load_data", usage="load IBD segment data into a redis database"
    )

    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file containing IBD segments. This file should have been generated with either GERMLINE or GERMLINE2",
    )

    parser.add_argument(
        "--min-cm",
        type=float,
        default=2.5,
        help="minimum centimorgan threshold to use in the ersa analysis. This value will be used to filter the ibd segments for the analysis",
    )

    parser.add_argument(
        "--port",
        type=int,
        default=6379,
        help="This is the port that will be used to connect to the redis database.",
    )

    args = parser.parse_args()
    # we need to create connect to the redis database
    #
    # Lets initialize the runtime state
    state = RuntimeState(args.min_cm)

    print(f"connecting to a redis database at localhost on port {args.port}")

    retries = Retry(ExponentialBackoff(), 8)
    try:
        r = redis.Redis(
            host="localhost",
            port=args.port,
            decode_responses=True,
            health_check_interval=3,
            retry=retries,
        )

        if r.ping():
            print("successfully connected to the Redis database")
        else:
            print("ERROR: the service is not available. Exiting program...")
            sys.exit(1)

    except redis.exceptions.ConnectionError as e:
        print(f"Failed to open the database becaause of the following error:\n{str(e)}")
        sys.exit(1)

    # Start reading in the file. The runtime state obj will have the centimoran threshold already
    file_iterator = read_segment_data(args.input, state)

    for segment_info in file_iterator:
        add_values_to_db(segment_info, r)

    # if segment_dict:  # Only if we have segments loaded
    #     additional_options = update_min_cm(
    #         additional_options, segment_dict
    #     )  # This will update min_cm passed into ersa if applicable


if __name__ == "__main__":
    main()
