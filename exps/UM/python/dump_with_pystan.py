import pystan
import pyreadr
from pyprojroot import here
from pathlib import PurePath


## cannot allocate memory error.
mydata = pyreadr.read_r(here(PurePath("data", "UM", "rstan",
                                      "sc_genewise.RData")).as_posix(),
                        use_objects=["N"])
