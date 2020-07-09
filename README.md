# ceRNA_multiMIR_starbase
---
title: "从mRNA到ceRNA network"
output: html_document

---



### 0.R包安装


```r
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/"))
options(download.file.method = 'libcurl')
options(url.method='libcurl')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require(multiMiR))BiocManager::install("multiMiR",ask = F,update = F)
library(multiMiR)
```

### 1.从mRNA得到miRNA

输入数据的获得方式：使用某个癌症的RNA-seq数据经过差异分析、PPI网络等的筛选，得到几个关键的mRNA，保存在9hubgenes.txt。


```r
x = read.table("9hubgenes.txt",stringsAsFactors = F)$V1;x
```

```
## [1] "MMP9"  "CXCL8" "ACTB"  "TGB1"  "STAT1" "TOP2A"
## [7] "CDK1"  "GNMT"  "ABAT"
```

我翻阅了很多网页工具和教程，发现multiMiR这个R包很优秀，结合了14个数据库，其中包括了有实验方法验证互作关系的mirTarbase。详见：![microRNAs靶基因数据库哪家强]https://mp.weixin.qq.com/s/n_UncYeGIQFLneTMK2rTXQ

这个工具可以从mRNA得到miRNA，也可以从miRNA得到mRNA，在这里使用前者，table = 'validated'这个参数是默认的，写上是为了和table = 'predicted'区分开，两者对应的数据库源不同，前者是经实验验证的，数量会更少一些，也更可靠一些。


```r
gene2mir <- get_multimir(org     = 'hsa',
                         target  = x,
                         table   = 'validated',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)
```

```
## Searching mirecords ...
## Searching mirtarbase ...
## Searching tarbase ...
```

```r
mit = gene2mir@data[gene2mir@data$database=="mirtarbase",];dim(mit)
```

```
## [1] 397  10
```

这里大材小用一下，只选了3个数据库中的一个，实验验证也分几个等级,最严谨的是Luciferase reporter assay，只选出它：

```r
table(mit$support_type)
```

```
## 
##        Functional MTI Functional MTI (Weak) 
##                    43                   353 
##    Non-Functional MTI 
##                     1
```

```r
mit = mit[stringr::str_detect(mit$experiment,
                            "Luciferase reporter assay"),];dim(mit)
```

```
## [1] 25 10
```

```r
miRNAs = unique(mit$mature_mirna_id)
```


### 2.从miRNA得到lncRNA

我查了一下相关的文献，lncRNA 和miRNA互作的数据库使用比较多的有三个：mircode，starbase，mirnet，我都看了一下，最后感觉starbase表现最好。

我在网页上戳戳戳了半天，发现starbase只能一次搜索一个miRNA，我疯了。想要下载它的lncRNA - miRNA 互作数据自己探索，没有在网页上找到直接下载的按钮，但找到了关于API的说明里有：
![](https://upload-images.jianshu.io/upload_images/9475888-f28a988ffff0597e.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
截图出自：http://starbase.sysu.edu.cn/tutorialAPI.php

在linux命令行用curl，写对筛选要求就可以获取数据啦。全部lncRNA - miRNA 互作数据的获取代码是：


```r
curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=lncRNA&miRNA=all&clipExpNum=0&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=all&cellType=all' > starBaseV3_hg19_CLIP-seq_lncRNA_all.txt &
```

下载的这句代码来自：https://www.jianshu.com/p/b7e4830c0b01

有了这个数据，读进R语言就可以随便玩耍啦。

```r
starbase = data.table::fread("starBaseV3_hg19_CLIP-seq_lncRNA_all.txt");dim(starbase)
```

```
## [1] 71952    17
```

很多文献里会把GENECODE里没有注释的lncRNA去掉，这个也很好操作，anno.Rdata来自于genecodev22版本的gtf文件，参考：![从TCGA表达矩阵中分别提取mRNA和lncRNA](https://mp.weixin.qq.com/s/bGoUbLuBdPteo-oG8ckMVw)

在今天的资料里这个anno.Rdata文件直接提供，可以反复使用的。


```r
load("anno.Rdata")
lnc_anno$gene_id = stringr::str_remove(lnc_anno$gene_id,"\\.\\d")
p1 = starbase$geneName %in% lnc_anno$gene_name;table(p1)
```

```
## p1
## FALSE  TRUE 
## 43924 28028
```

```r
p2 = starbase$geneID %in% lnc_anno$gene_id;table(p2)
```

```
## p2
## FALSE  TRUE 
## 11016 60936
```

```r
starbase = starbase[p2,];dim(starbase)
```

```
## [1] 60936    17
```

```r
lnc_mi = starbase[starbase$miRNAname %in% miRNAs,]
colnames(lnc_mi)
```

```
##  [1] "miRNAid"      "miRNAname"    "geneID"      
##  [4] "geneName"     "geneType"     "chromosome"  
##  [7] "start"        "end"          "strand"      
## [10] "clipExpNum"   "degraExpNum"  "RBP"         
## [13] "merClass"     "miRseq"       "align"       
## [16] "targetSeq"    "pancancerNum"
```

对starbase数据库提供的数据列名的解释里，有两个比较重要的：

> **clipExpNum** ：The number of CLIP-seq experiments ; **pancancerNum** : Number of Cancer types (Pan-Cancer) that miRNA-target show anti-correlation relationships (pearson  correlation: r<0, p-value<0.05).

所以pancerNum是几 就意味着这对miRNA-lncRNA在多少种癌症中表达量负相关，省掉了很多计算。

数据库的tutorial里面还提到了一个比较严格的筛选标准：
CLIP evidence (>=5), degradome evidence (>=1), Pan-Cancer (>=10), program number (>=5) and predicted program (None).

degradome evidence 的限制条件，加上和不加数量有差别。


```r
p2 = lnc_mi$pancancerNum >10  &lnc_mi$clipExpNum>4;table(p2)
```

```
## p2
## FALSE  TRUE 
##  2143    19
```

```r
p3 = lnc_mi$pancancerNum >10  &lnc_mi$clipExpNum>4 & lnc_mi$degraExpNum >0;table(p3)
```

```
## p3
## FALSE  TRUE 
##  2160     2
```

```r
lnc_mi$geneName[p2]
```

```
##  [1] "NORAD"      "MALAT1"     "MALAT1"    
##  [4] "SNHG1"      "AC105339.2" "KCNQ1OT1"  
##  [7] "NORAD"      "AC090198.1" "NORAD"     
## [10] "NORAD"      "NORAD"      "GAS5"      
## [13] "FGD5-AS1"   "PITPNA-AS1" "AC078846.1"
## [16] "OIP5-AS1"   "TUG1"       "TUG1"      
## [19] "HAGLR"
```

```r
lnc_mi$geneName[p3]
```

```
## [1] "SNHG1"      "PITPNA-AS1"
```

这里就不加上degradome evidence 的限制条件了。

### 3.准备cytoscape的输入文件


```r
ez1 = mit[,3:4]
ez2 = lnc_mi[p2,c(2,4)]
colnames(ez1) = colnames(ez2)

library(dplyr)
ez = rbind(ez1,ez2);dim(ez)
```

```
## [1] 44  2
```

```r
ez = distinct(ez,miRNAname,geneName);dim(ez)
```

```
## [1] 38  2
```

网络的颜色与形状需要按照RNA种类来区分，所以制作一个RNA与种类对应关系的数据框。


```r
tp = data.frame(nodes = c(ez1$miRNAname,
                          ez2$miRNAname,
                          ez1$geneName,
                          ez2$geneName),
                type = rep(c("mi","pc","lnc"),
                           times = c(nrow(ez1)+nrow(ez2),
                                     nrow(ez1),
                                     nrow(ez2))
                           )
                )
dim(tp);head(tp)
```

```
## [1] 88  2
```

```
##             nodes type
## 1    hsa-miR-644a   mi
## 2   hsa-miR-31-5p   mi
## 3  hsa-miR-23a-3p   mi
## 4   hsa-miR-93-5p   mi
## 5 hsa-miR-302c-3p   mi
## 6 hsa-miR-302d-3p   mi
```

```r
tp = distinct(tp,nodes,.keep_all = T)
table(tp$type)
```

```
## 
## lnc  mi  pc 
##  13  23   5
```

最后是导出。


```r
write.table(ez,file = "ez.txt",quote = F,sep = "\t",row.names = F)
write.table(tp,file = "tp.txt",quote = F,sep = "\t",row.names = F)
```

