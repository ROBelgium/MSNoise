import os
import pydoc
import pandas as pd
from obspy.core.util import AttribDict

def get_defaults():
    df = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), "default.csv"),
                     delimiter=",", encoding="latin-1", index_col=0)
    df["type"] = [pydoc.locate(t) for t in df["type"].values]
    df = df.fillna("")

    default = AttribDict()
    for id, row in df.iterrows():
        if len(row.possible_values):
            df.loc[id].definition += " " + row.possible_values.replace(row.default, "[%s]"%row.default)
        # elif len(row.default):
        #     df.loc[id].definition += " " + "[%s]" % row.default
        default[id] = AttribDict(row.to_dict())
    return default

default = get_defaults()
default_datetime_fields = ["startdate", "enddate", "ref_begin", "ref_end"]
