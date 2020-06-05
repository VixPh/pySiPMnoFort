from pyspark.sql import SparkSession
from main import *

spark = SparkSession.builder.appName('pySiPM').getOrCreate()
sc = spark.sparkContext

sc.addFile('files.zip')
