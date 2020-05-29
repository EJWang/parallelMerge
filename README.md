# paralleMerge

> A tool help to merge RNA data based on metadata JSON

> calculate the number of normal, tumor and generate the merged RNA Matrix text file

> put it on your data folder with metadata then run it, no need to decompress

## Usage
Usage: python3 parallelMerge.py <metadata JSON>

## Example of shell output

```shell
Current: c4a29253-97ba-40d8-b920-1b90871d94b3/49025bb6-3da7-4e27-b8a3-c3815ede80eb.FPKM.txt.gz......
Current: 4303d5d3-3b74-47f9-8b39-6bc1493c2661/571143cf-c9f4-4490-8ff7-084ab5b297d3.FPKM.txt.gz......

---------------- Stats ----------------
Normal Count: 104
Tumor Count: 828
------------------------------------------

Result has been generated to mRNA_Matrix.txt

Process done, it takes: 19.75sec


Please wait, program exiting.......
```
