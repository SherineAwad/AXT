import pyreadr
import pandas as pd

result = pyreadr.read_r("model_mouse_full.rda")
df = result['model_mouse_full']

# Save as CSV
df.to_csv("progeny_mouse.csv", index=False)

print(df.head())
print(f"\nShape: {df.shape}")
