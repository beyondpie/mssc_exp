import pystan
import pyreadr
import argparse
from typing import Dict
import pandas as pd
import copy
import logging
from pyprojroot import here
from pathlib import PurePath
import numpy as np

logging.basicConfig(format='%(asctime)s %(message)s')


def to_numpy(mydata: Dict) -> Dict:
    tmp = copy.deepcopy(mydata)
    for name, value in tmp.items():
        if isinstance(value, (pd.Series, pd.DataFrame)):
            v = value.to_numpy()
            if v.shape == (1,1):
                tmp[name] = int(v[0,0])
                # if v.dtype == np.float:
                #     logging.warnings(f"{name} shape is (1,1) float, change it to float")
                #     tmp[name] = float(v[0,0])
                # elif v.dtype == np.int:
                #     logging.warnings(f"{name} shape is (1,1) int, change it to int")
                #     tmp[name] = int(v[0,0])
                # else:
                #     logging.error(f"{name} shape is (1,1) {v[0,0]} with unknown type")
            else:
                tmp[name] = v
    return tmp


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "transform RData for stan using pystan rdump.")
    parser.add_argument("-d",
                        "--data_dir",
                        type=str,
                        default="data",
                        help="data dir")
    parser.add_argument("-s",
                        "--sub",
                        type=str,
                        default="UM",
                        help="sub dir under data")
    parser.add_argument("-i", "--infl", type=str, help="input file")
    parser.add_argument("-o", "--outfl", type=str, help="ouput file")
    args = parser.parse_args()

    inputf = here(PurePath(args.data_dir, args.sub, args.infl)).as_posix()
    mydata: Dict = pyreadr.read_r(inputf)
    npdata: Dict = to_numpy(mydata)

    outputf = here(PurePath(args.data_dir, args.sub, args.outfl)).as_posix()
    pystan.misc.stan_rdump(npdata, outputf)
    logging.info("rdump done.")
