import geopandas as geopandas
import pandas as pd
import geopandas as gpd
from libpysal.weights import Queen
import numpy as np
from esda.moran import Moran

name=["UKF241_C", "UKF242_C", "UKF242_T", "UKF243_T", "UKF248_C", "UKF248_T", "UKF251_T", "UKF255_T", "UKF256_C", "UKF259_C", "UKF259_T", "UKF260_T",
    "UKF262_T", "UKF265_T", "UKF266_T", "UKF268_T", "UKF269_T", "UKF270_T", "UKF275_T", "UKF296_T", "UKF304_T", "UKF313_T", "UKF334_C", "UKF334_T"]
index1 = [];index2 = [];index3=[]

for i in range(0,24):
    # x-coordinates of [0,100]×[0,100] meshgrids flattened in column-major order
    X = np.loadtxt(open( str(name[i]) + "/position.csv",
        "rb"), delimiter=",", skiprows=1, usecols=[0])
    # y-coordinates of [0,100]×[0,100] meshgrids flattened in column-major order
    Y = np.loadtxt(open( str(name[i]) + "/position.csv",
        "rb"), delimiter=",", skiprows=1, usecols=[1])
    Z1 = np.loadtxt(
        open(str(name[i]) + "/malignant.csv", "rb"),
        delimiter=",", skiprows=1, usecols=[0])   # malignant
    Z2 = np.loadtxt(
        open( str(name[i]) + "/endothelial.csv", "rb"),
        delimiter=",", skiprows=1, usecols=[0])  # endothelial
    Z3 = np.loadtxt(open( str(name[i]) + "/phenotype.csv",
        "rb"), delimiter=",", skiprows=1, usecols=[0])  # phenotype

    dict = {
         'x': X,
         'y': Y,
         'z': Z1}
    df = pd.DataFrame(dict)
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.x, df.y))

    w = Queen.from_dataframe(gdf)
    y =gdf.z
    moran= Moran(y,w)
    index1.append(moran.I)

    dict = {
        'x': X,
        'y': Y,
        'z': Z2}
    df = pd.DataFrame(dict)
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.x, df.y))

    w = Queen.from_dataframe(gdf)
    y = gdf.z
    moran = Moran(y, w)
    index2.append(moran.I)

    dict = {
        'x': X3,
        'y': Y3,
        'z': Z3}
    df = pd.DataFrame(dict)
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.x, df.y))

    w = Queen.from_dataframe(gdf)
    y = gdf.z
    moran = Moran(y, w)
    index3.append(moran.I)

test1=pd.DataFrame(
    {
        'MoranI(m)':index1,
        'MoranI(e)':index2,
        'MoranI(p)':index3
    }
)
test1.to_csv('MoranIndex.csv')



