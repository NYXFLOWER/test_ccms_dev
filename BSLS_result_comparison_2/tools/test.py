import pandas as pd
import numpy as np
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

df = pd.read_csv(input_path, sep="\t")
df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
df["c"] = [7, 8, 9]
df.to_csv(output_path, index=False)
