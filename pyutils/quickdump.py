import pystan
import pyreadr
import argparse
from typing import Dict
from numbers import Number
from numpy import np
import pandas as pd
import copy
import logging

logging.basicConfig(format='%(asctime)s %(message)s')


def test():
    from pyprojroot import here
    from pathlib import PurePath

    mydata = pyreadr.read_r(here(
        PurePath("data", "UM", "rstan", "sc_genewise.RData")).as_posix(),
                            use_objects=["N"])
    mydata = pyreadr.read_r(
        here(PurePath("data", "UM", "rstan", "test.RData")).as_posix())
    Xcg = {'Xcg': mydata['Xcg'].to_numpy()}
    pystan.misc.stan_rdump(
        Xcg, here(PurePath("data", "UM", "rstan", "test.pydump").as_posix()))


def to_numpy(mydata: Dict) -> Dict:
    tmp = copy.deepcopy(mydata)
    for name, value in tmp.items():
        if isinstance(value, (pd.Series, pd.DataFrame)):
            logging.info(f"{name} is pandas format, now change it to numpy")
            tmp[name] = value.to_numpy()
    return tmp


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "transform RData for stan using pystan rdump.")
    parser.add_argument("-i", "--infl", type=str, help="input file")
    parser.add_argument("-o", "--outfl", type=str, help="ouput file")
    args = parser.parse_args()
    mydata: Dict = pyreadr.read_r(args.infl)
    logging.info("detect pandas format, and transform to numpy ...")
    npdata: Dict = to_numpy(mydata)
    logging.info("passing to pystan rdump...")
    pystan.misc.stan_rdump(npdata, args.outfl)
    logging.info("rdump done.")
