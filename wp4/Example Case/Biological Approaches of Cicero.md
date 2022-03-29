# Biological approaches of Cicero

## Working with Cicero co-accessibility-scores (Research Question 1)


The main output of `Cicero` is the pairwise co-accessibility score. This information is stored in the file ***conns.csv***. Within it lays a data frame with two peaks and one assigned coaccess score. This score ranges from -1 to 1. All peaks originate from a single tissue sample, in this case ***esophagus_mucosa***. 

This work package wants to look at the scores in depth and both sort and assign the peaks. `Cicero` is able to assign peaks to known biological structures, as promoter regions. The information that links peaks with promoter regions of genes is provided in the file ***cds_sites.csv***. 

Once the peaks are assigned to a gene it would be very interesting to compare the peaks of one gene. Which leads to the first research question:
### Exists an open promoter region in two cell types with different accessibility of distal elements?
This section will get to the bottom of this question.

In order to find distal peaks to peaks assigned to promoters, one has to connect the two data frames ***conns.csv*** and ***cds_sites.csv*** with one antoher. This ensures that only peaks with promoters assigned come into account. 
To create this new data frame the files have to be loaded first. 


```R
library(dplyr)
library(stringr)
library(Matrix)
```


***conns.csv*** is loaded to the variable ***conns***, by the function `read.csv()`, specifing which seperator was used. The columnames should be ***Peak1/2*** and ***coaccess***. 


```R
conns <- read.csv("/mnt/workspace_stud/stud10/output/esophagus_mucosa/conns.csv", 
                      sep = ',', row.names = NULL)
conns <- conns[c("Peak1", "Peak2", "coaccess")]
tail(conns)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>11837165</th><td>chr7_99990880_99991280</td><td>chr7_100243424_100243824</td><td> 2.418531e-04</td></tr>
	<tr><th scope=row>11837166</th><td>chr7_99990880_99991280</td><td>chr7_100247994_100248394</td><td>-3.751595e-05</td></tr>
	<tr><th scope=row>11837167</th><td>chr7_99990880_99991280</td><td>chr7_100180677_100181077</td><td> 3.557228e-04</td></tr>
	<tr><th scope=row>11837168</th><td>chr7_99990880_99991280</td><td>chr7_100219198_100219598</td><td> 2.848364e-04</td></tr>
	<tr><th scope=row>11837169</th><td>chr7_99990880_99991280</td><td>chr7_100222053_100222453</td><td> 4.394876e-05</td></tr>
	<tr><th scope=row>11837170</th><td>chr7_99990880_99991280</td><td>chr7_100240015_100240415</td><td>-2.598716e-05</td></tr>
</tbody>
</table>



***cds_site.csv*** is loaded to the variable ***cds_site***. To connect this two data frames with one another one column has to have the same name, in this case ***Peak2***. 


```R
cds_site <- read.csv("/mnt/workspace_stud/stud10/output/esophagus_mucosa/cds_sites.csv", 
                      sep = ',', row.names = NULL)
colnames(cds_site) <- c("num", "Peak2", "gene")
cds_site <- cds_site[,c("Peak2", "gene")]
head(cds_site)
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Peak2</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr1_817141_817541</td><td>FAM87B </td></tr>
	<tr><th scope=row>2</th><td>chr1_925472_925872</td><td>SAMD11 </td></tr>
	<tr><th scope=row>3</th><td>chr1_940684_941084</td><td>SAMD11 </td></tr>
	<tr><th scope=row>4</th><td>chr1_970716_971116</td><td>PLEKHN1</td></tr>
	<tr><th scope=row>5</th><td>chr1_973139_973539</td><td>PLEKHN1</td></tr>
	<tr><th scope=row>6</th><td>chr1_982026_982426</td><td>PERM1  </td></tr>
</tbody>
</table>



Create new data frame ***dff*** by comparison of ***conns.Peak2*** and ***cds_site.Peak2***. Afterwards only peaks found in both ***Peak2*** columns are in the final data frame. ***dff*** now displays ***Peak1***, ***Peak2*** and their calculated ***co-access***. ***Peak2*** is also asigned to a ***gene***.  


```R
# combine co-access table with cds site table to assign promoters to peaks
dff <- inner_join(conns, cds_site, by = "Peak2")
# remove entries with NA
dff <- dff[complete.cases(dff),]
# order data by gene with inner order by coaccess
dff <- dff[order(dff$gene, dff$coaccess),]
head(dff)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>138115</th><td>chr17_35766915_35767315</td><td>chr17_35889180_35889580</td><td>-0.03748543</td><td>7SK</td></tr>
	<tr><th scope=row>138139</th><td>chr17_35774105_35774505</td><td>chr17_35889180_35889580</td><td>-0.03729434</td><td>7SK</td></tr>
	<tr><th scope=row>138055</th><td>chr17_35748257_35748657</td><td>chr17_35889180_35889580</td><td>-0.03530089</td><td>7SK</td></tr>
	<tr><th scope=row>138183</th><td>chr17_35778095_35778495</td><td>chr17_35889180_35889580</td><td>-0.03513867</td><td>7SK</td></tr>
	<tr><th scope=row>138047</th><td>chr17_35741024_35741424</td><td>chr17_35889180_35889580</td><td>-0.02839489</td><td>7SK</td></tr>
	<tr><th scope=row>138219</th><td>chr17_35793196_35793596</td><td>chr17_35889180_35889580</td><td>-0.02617012</td><td>7SK</td></tr>
</tbody>
</table>



The ***dff*** data frame has to be filtered into a subset data frame, only using scores higher than 0.2 and smaller than -0.01. 


```R
dff <- subset(dff, coaccess > 0.2 | coaccess < -0.01)
head(dff)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>138115</th><td>chr17_35766915_35767315</td><td>chr17_35889180_35889580</td><td>-0.03748543</td><td>7SK</td></tr>
	<tr><th scope=row>138139</th><td>chr17_35774105_35774505</td><td>chr17_35889180_35889580</td><td>-0.03729434</td><td>7SK</td></tr>
	<tr><th scope=row>138055</th><td>chr17_35748257_35748657</td><td>chr17_35889180_35889580</td><td>-0.03530089</td><td>7SK</td></tr>
	<tr><th scope=row>138183</th><td>chr17_35778095_35778495</td><td>chr17_35889180_35889580</td><td>-0.03513867</td><td>7SK</td></tr>
	<tr><th scope=row>138047</th><td>chr17_35741024_35741424</td><td>chr17_35889180_35889580</td><td>-0.02839489</td><td>7SK</td></tr>
	<tr><th scope=row>138219</th><td>chr17_35793196_35793596</td><td>chr17_35889180_35889580</td><td>-0.02617012</td><td>7SK</td></tr>
</tbody>
</table>



The variable ***oc_dff*** stores the joined ***conns*** data frame with the difference that the values in the column ***coaccess*** are replaced by '-' if the value is negative and '+' if the value is positiv. The values itself do not matter at this point, only if the peak is open with the promoter region or closed. 


```R
oc_dff <- dff
oc_dff[oc_dff < 0] <- '-'
oc_dff$coaccess[oc_dff$coaccess > 0] <- '+'
head(oc_dff)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>138115</th><td>chr17_35766915_35767315</td><td>chr17_35889180_35889580</td><td>-</td><td>7SK</td></tr>
	<tr><th scope=row>138139</th><td>chr17_35774105_35774505</td><td>chr17_35889180_35889580</td><td>-</td><td>7SK</td></tr>
	<tr><th scope=row>138055</th><td>chr17_35748257_35748657</td><td>chr17_35889180_35889580</td><td>-</td><td>7SK</td></tr>
	<tr><th scope=row>138183</th><td>chr17_35778095_35778495</td><td>chr17_35889180_35889580</td><td>-</td><td>7SK</td></tr>
	<tr><th scope=row>138047</th><td>chr17_35741024_35741424</td><td>chr17_35889180_35889580</td><td>-</td><td>7SK</td></tr>
	<tr><th scope=row>138219</th><td>chr17_35793196_35793596</td><td>chr17_35889180_35889580</td><td>-</td><td>7SK</td></tr>
</tbody>
</table>



Select yourself a gene you want to have a closer look on. You can either choose from the list provided by the variable ***list_genes*** or you can look through the data frame ***oc_dff***. ***list_genes*** is a sorted data frame of the genes and its occurrences in the modified coaccess ***oc_dff*** data frame. 

The genes are assigned to the peaks of ***Peak2***, the ***Freq*** column of ***list_genes*** describes then the occurrences of a gene to be assigned to a peak. 


```R
list_genes <-sort(table(unlist(oc_dff$gene)))
list_genes <- as.data.frame(list_genes)
colnames(list_genes) <- c("gene", "Freq")
list_genes <- list_genes[list_genes$Freq>30, ]
tail(list_genes)
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>gene</th><th scope=col>Freq</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>2326</th><td>CD27   </td><td>161</td></tr>
	<tr><th scope=row>2327</th><td>TIGIT  </td><td>174</td></tr>
	<tr><th scope=row>2328</th><td>CD6    </td><td>177</td></tr>
	<tr><th scope=row>2329</th><td>SLC1A2 </td><td>179</td></tr>
	<tr><th scope=row>2330</th><td>ST6GAL1</td><td>193</td></tr>
	<tr><th scope=row>2331</th><td>SEPTIN9</td><td>283</td></tr>
</tbody>
</table>



Create a sublist of ***oc_dff*** with all rows, which gene value is (as random example) *ST6GAL1*. 


```R
gene_dff <- oc_dff[which(oc_dff$gene == 'ST6GAL1'), ]
head(gene_dff)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>449001</th><td>chr3_186949177_186949577</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td></tr>
	<tr><th scope=row>449421</th><td>chr3_187173738_187174138</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td></tr>
	<tr><th scope=row>448822</th><td>chr3_186758132_186758532</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td></tr>
	<tr><th scope=row>449538</th><td>chr3_187352677_187353077</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td></tr>
	<tr><th scope=row>449471</th><td>chr3_187207599_187207999</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td></tr>
	<tr><th scope=row>448865</th><td>chr3_186808718_186809118</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td></tr>
</tbody>
</table>



Read in the data frame which combines ***Peak2*** with the clusters. This data frame will be stored in the variable ***cluster_peaks***. 
This data frame contains the peaks that are similar to the values of ***Peak2*** of ***gene_dff***. And the counts of the clusters. In this example data of the esophagus_mucosa there were 17 clusters assigned. For the gene *ST6CAL1* only 5 different Peaks were assigned. 


```R
cluster_peaks <- read.csv("/mnt/workspace_stud/stud10/cluster_peaks.csv", 
                      sep = ',', row.names = NULL)
cluster_peaks <- cluster_peaks[complete.cases(cluster_peaks), ]
colnames(cluster_peaks) <- c("Peak2", "0", "1", "2", "3", "4", "5", 
                             "6", "7", "8", "9", "10", "11", "12",
                            "13", "14", "15", "16", "17")
head(cluster_peaks)
```

<table class="dataframe">
<caption>A data.frame: 5 × 19</caption>
<thead>
	<tr><th></th><th scope=col>Peak2</th><th scope=col>0</th><th scope=col>1</th><th scope=col>2</th><th scope=col>3</th><th scope=col>4</th><th scope=col>5</th><th scope=col>6</th><th scope=col>7</th><th scope=col>8</th><th scope=col>9</th><th scope=col>10</th><th scope=col>11</th><th scope=col>12</th><th scope=col>13</th><th scope=col>14</th><th scope=col>15</th><th scope=col>16</th><th scope=col>17</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr3_186985931_186986331</td><td>2</td><td>1</td><td>21</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>1</td><td>22</td><td>0</td><td>1</td><td>1</td><td> 9</td><td>3</td><td>0</td><td>2</td><td>2</td></tr>
	<tr><th scope=row>2</th><td>chr3_186996260_186996660</td><td>2</td><td>1</td><td> 6</td><td>2</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>40</td><td>1</td><td>0</td><td>1</td><td>11</td><td>0</td><td>0</td><td>2</td><td>2</td></tr>
	<tr><th scope=row>3</th><td>chr3_187021630_187022030</td><td>0</td><td>0</td><td> 1</td><td>3</td><td>1</td><td>0</td><td>1</td><td>1</td><td>2</td><td>13</td><td>1</td><td>0</td><td>0</td><td>20</td><td>0</td><td>0</td><td>3</td><td>0</td></tr>
	<tr><th scope=row>4</th><td>chr3_187024500_187024900</td><td>2</td><td>1</td><td>16</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td> 7</td><td>0</td><td>0</td><td>0</td><td> 2</td><td>1</td><td>0</td><td>2</td><td>2</td></tr>
	<tr><th scope=row>5</th><td>chr3_187025183_187025583</td><td>0</td><td>1</td><td> 4</td><td>2</td><td>0</td><td>2</td><td>0</td><td>0</td><td>3</td><td> 8</td><td>0</td><td>0</td><td>0</td><td>21</td><td>1</td><td>0</td><td>4</td><td>0</td></tr>
</tbody>
</table>



Now the column ***max*** is generated by hand. This gives information about the cluster/s that should be assigned to a peak. 
And the data frame is reduced to the two columns ***Peak2*** and ***max***


```R
cluster_peaks$max <- c("9_2", "9", "13_20", "2", "13")
cluster_peaks <- cluster_peaks[, c("Peak2", "max")]
```

`inner_join()` combines now the data of ***gene_dff*** and ***cluster_peaks*** by the column ***Peak2***. The combined data frame ***cluster_gene_dff*** is sorted by ***Peak2*** and within this assortment by ***coaccess***. 


```R
cluster_gene_dff <- inner_join(gene_dff, cluster_peaks, by = "Peak2")
cluster_gene_dff <- cluster_gene_dff[order(cluster_gene_dff$Peak2, 
                                           cluster_gene_dff$coaccess),]
head(cluster_gene_dff)
```


<table class="dataframe">
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th><th scope=col>max</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>23</th><td>chr3_186949177_186949577</td><td>chr3_186985931_186986331</td><td>-</td><td>ST6GAL1</td><td>9_2</td></tr>
	<tr><th scope=row>32</th><td>chr3_186640920_186641320</td><td>chr3_186985931_186986331</td><td>-</td><td>ST6GAL1</td><td>9_2</td></tr>
	<tr><th scope=row>33</th><td>chr3_187207599_187207999</td><td>chr3_186985931_186986331</td><td>-</td><td>ST6GAL1</td><td>9_2</td></tr>
	<tr><th scope=row>34</th><td>chr3_187173738_187174138</td><td>chr3_186985931_186986331</td><td>-</td><td>ST6GAL1</td><td>9_2</td></tr>
	<tr><th scope=row>36</th><td>chr3_186758132_186758532</td><td>chr3_186985931_186986331</td><td>-</td><td>ST6GAL1</td><td>9_2</td></tr>
	<tr><th scope=row>41</th><td>chr3_187174229_187174629</td><td>chr3_186985931_186986331</td><td>-</td><td>ST6GAL1</td><td>9_2</td></tr>
</tbody>
</table>



Sublist ***cluster_gene_dff*** by cluster, by checking the unique values in the clusters, this leads to a list of the different clusters in the data frame ***cluster_gene_dff***. This list is stored in the variable ***list_max_values***. The length of this list (***len***) is the number of different clusters, in this case, the number of different subset that should be created. 


```R
list_max_values <- unique(cluster_gene_dff$max)
len <- length(list_max_values)
```

This code chunk generates a list ***eq*** with the lists ***lhs*** and ***rhs*** put together in a string. This string then initializes the variables of the different sublists (names created in ***lhs***) to the entries of ***list_max_value*** (clusters, created in ***rhs***). 


```R
lhs  <- paste("sub",    1:len,     sep="")
rhs  <- paste("list_max_values[[",1:len,"]]", sep="")
eq   <- paste(paste(lhs, rhs, sep="<-"), collapse=";")
eval(parse(text=eq))
```

Now the assignments are stored in the variables ***subX***. This code then assigns the whole data frame matching the condition in variable ***subX*** to the variables ***sublistX***. 


```R
lhs  <- paste("sublist",    1:len,     sep="")
rhs  <- paste("cluster_gene_dff[which(cluster_gene_dff$max == sub",1:len, "),]", sep="")
eq   <- paste(paste(lhs, rhs, sep="<-"), collapse=";")
eval(parse(text=eq))
# create list with sublists 
sublist <- list(sublist1, sublist2, sublist3, sublist4, sublist5)
```

*Check 1: Are there same Peaks with differing co-access allocation*: 

This iteration goes though the list of the sublists and compares the ***ocaccess*** columns with one another. And if it finds a Peak that is open in one cluster, but not in the other it will print the number of not matching sublists, or will tell that the data frame which was compared is equal with every other frame. 


```R
counter1 <- 0
counter2 <- 0
for (a in sublist){
    counter1 <<- counter1 + 1
    for (b in sublist){
        counter2 <<- counter2 + 1
        ij_list <- inner_join(a, b, by = "Peak1")
        ij_list <- ij_list[, c("Peak1", "coaccess.x", "max.x", "coaccess.y", "max.y")]
        x <- ij_list$coaccess.x
        y <- ij_list$coaccess.y
        if (! identical(x, y)){
            print(paste(counter1, "not equal to ", counter2))
        }
        else {
        }
    
    }
    counter2 <- 0
    print(paste(counter1, " is equal with every data frame."))
}
```

    [1] "1  is equal with every data frame."
    [1] "2  is equal with every data frame."
    [1] "3  is equal with every data frame."
    [1] "4  is equal with every data frame."
    [1] "5  is equal with every data frame."

So no Peaks were found with differing co-accessibility sheme. 

*2. Check: Are there unique Peaks with co-accessibility scores found per Cluster:* 

This question can be answered by comparing all sublists and creating a data frame with just uniquely found peaks in the column ***peak1*** that are just found in one cluster. 
***freq*** stores the a table of occurrences of values in ***Peak1***. This tables is sublisted to only take peaks with occurrences of 1. The colnames have to be adjusted in order to conncect ***freq*** with ***cluster_gene_dff***. 


```R
freq <- as.data.frame(table(cluster_gene_dff$Peak1))
freq <- freq[which(freq$Freq == 1), ]
colnames(freq) <- c("Peak1", "Freq")
head(freq)
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Freq</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr3_186521512_186521912</td><td>1</td></tr>
	<tr><th scope=row>18</th><td>chr3_186882812_186883212</td><td>1</td></tr>
	<tr><th scope=row>19</th><td>chr3_186909587_186909987</td><td>1</td></tr>
	<tr><th scope=row>25</th><td>chr3_186952986_186953386</td><td>1</td></tr>
	<tr><th scope=row>26</th><td>chr3_186971620_186972020</td><td>1</td></tr>
	<tr><th scope=row>36</th><td>chr3_187005831_187006231</td><td>1</td></tr>
</tbody>
</table>



`inner_join()` connects the data frames ***cluster_gene_dff*** and ***freq*** to a data frame called ***cluster_gene_dff_unique*** which only contains peaks uniquely found in one cluster. 


```R
cluster_gene_dff_unique <- inner_join(freq, cluster_gene_dff, by = "Peak1")
cluster_gene_dff_unique <- cluster_gene_dff_unique[order(cluster_gene_dff_unique$max, 
                                           cluster_gene_dff_unique$coaccess),]
cluster_gene_dff_unique
```


<table class="dataframe">
<caption>A data.frame: 13 × 6</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Freq</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th><th scope=col>max</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>5</th><td>chr3_186971620_186972020</td><td>1</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>10</th><td>chr3_187209513_187209913</td><td>1</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>11</th><td>chr3_187291526_187291926</td><td>1</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>12</th><td>chr3_187319222_187319622</td><td>1</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>13</th><td>chr3_187359798_187360198</td><td>1</td><td>chr3_187025183_187025583</td><td>-</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>3</th><td>chr3_186909587_186909987</td><td>1</td><td>chr3_187025183_187025583</td><td>+</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>4</th><td>chr3_186952986_186953386</td><td>1</td><td>chr3_187025183_187025583</td><td>+</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>6</th><td>chr3_187005831_187006231</td><td>1</td><td>chr3_187025183_187025583</td><td>+</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>9</th><td>chr3_187041339_187041739</td><td>1</td><td>chr3_187025183_187025583</td><td>+</td><td>ST6GAL1</td><td>13   </td></tr>
	<tr><th scope=row>2</th><td>chr3_186882812_186883212</td><td>1</td><td>chr3_187021630_187022030</td><td>-</td><td>ST6GAL1</td><td>13_20</td></tr>
	<tr><th scope=row>8</th><td>chr3_187025183_187025583</td><td>1</td><td>chr3_187021630_187022030</td><td>+</td><td>ST6GAL1</td><td>13_20</td></tr>
	<tr><th scope=row>1</th><td>chr3_186521512_186521912</td><td>1</td><td>chr3_186996260_186996660</td><td>+</td><td>ST6GAL1</td><td>9    </td></tr>
	<tr><th scope=row>7</th><td>chr3_187019571_187019971</td><td>1</td><td>chr3_186996260_186996660</td><td>+</td><td>ST6GAL1</td><td>9    </td></tr>
</tbody>
</table>



The data frame ***cluster_gene_dff_unique*** stores pairwise Peak co accessibility by clusters. The column ***Peak2*** is assigned to a promoter region (***gene***). ***Peak1*** is a peak uniquely found in the cluster the column ***max*** provides. The peaks in ***Peak1*** are distal to ***Peak2***. 

Especially the open Peaks are interesting: ***coaccess*** with '+'. The findings of these peaks could be a hint for gene expression within a cluster. What must be considered, however, is the fact that Cicero did not create any other coaccessible scores for these peaks. 

This could mean several things: 
1. Open promoter regions in other clusters do not provide gene expression because no distal peaks are coaccessbible for transcription factors or enhancers to bind.
2. open promoter regions in other peaks have other distal sites to which expression regulatory elements can attach
3. these distal peaks found are not related to expression but describe, for example, the coaccessibility of one gene to another.

***
*Optinal Code Chunk*: 
If cl_indata should be filtered before connecting clusters to peaks, this code helps the pre-filtering. 
A data frame is derived from ***dff*** with column ***peak2*** as rownames. ***cl_indata*** is now filtered by matching rownames with the derived data frame ***coacc_dff*** and stored in the directory of the tissue. 


```R
# create data frame with requirements to filter cl_indata 
coacc_dff <- gene_dff
rownames(coacc_dff) <- NULL
coacc_dff <- coacc_dff[!duplicated(coacc_dff$Peak2),]
row.names(coacc_dff) <- coacc_dff$Peak2
```


```R
# store filtered matrix in directory
write(x = rownames(subset), 
      file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/peaks.tsv")
write(x = colnames(subset), 
      file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cluster.tsv")

subset <- Matrix(subset , sparse = T )
writeMM(obj = subset, 
        file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cl_indata_esophagus_mucosa.mtx")
```

***
## Extension of Research Question 1: Working with gene activity scores

Accessibility at promoters is found to be a poor prediction of gene expression. Cicero links can be used to get a sense of an open promoter region and it's assigned distal sites. The combined score is called the "gene activity score". 

The previously created matrix is stored in `output/<tissue>/cicero_genes_activity.mtx`. The corresponding row names and col names are also stored there. To read the matrix the library(Matrix) is needed. 


```R
act <- Matrix::readMM("/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_genes_activity.mtx")
rows <- read.csv('/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_gene_activities_rows.tsv', 
                 col.names = 'rows', sep = ',', row.names = NULL, header = FALSE)
cols <- read.csv('/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_gene_activities_cols.tsv',  
                 col.names = 'barcode', sep = ',', header = FALSE)
```

To be able to assign clusters to the barcodes at the end, the data frame storing clusters and barcodes must be loaded into the variable ***bar_cluster***. This data frame must be modified to perform `inner_join()` with the column data frame ***cols*** of ***act***. 


```R
bar_cluster <- read.csv('/mnt/workspace_stud/stud10/output/esophagus_mucosa/bar_cluster.tsv',  
                 col.names = c("num", 'barcode'), sep = ',', header = FALSE)
bar_cluster[c('barcode', 'cluster')] <- str_split_fixed(bar_cluster$barcode, '_', 2)
bar_cluster_df <- bar_cluster[-1,]
bar_cluster_df$num <- NULL
```

`inner_join()` connects the data frames ***cols*** and ***bar_cluster_df*** by the column ***barcode***. In the next step a column is added that connects barcodes and clusters with the character _. 


```R
barcode_dff <- inner_join(cols, bar_cluster_df , by = "barcode")
barcode_dff$joined <- paste(barcode_dff$barcode, barcode_dff$cluster, sep="_")
```

Row and col names of ***act*** can now be set. 


```R
# set col and rown names of act.
row.names(act) <- rows$rows
colnames(act) <- barcode_dff$joined
```

Subset gene activity score matrix that only the selected gene ***ST6GAL1*** is displayed. This new data frame is then connected with the barcodes and clusters, filtered to only values > 0.1 and ordered by the scores. 


```R
include_list <- c("ST6GAL1")
gene_matrix <- suppressWarnings(subset(act, rownames(act) %in% include_list))
gene_matrix <- cbind(gene_matrix, barcode_dff$joined)
colnames(gene_matrix) <- c("activity_score", "barcode_cluster")

gene_matrix <- as.data.frame(gene_matrix)
gene_matrix[c('barcode','cluster')] <- str_split_fixed(gene_matrix$barcode_cluster, '_', 2)
gene_matrix <- gene_matrix[gene_matrix$activity_score > 0.1, ]
gene_matrix <- gene_matrix[, c("activity_score", "cluster")]
gene_matrix <- gene_matrix[order(gene_matrix$activity_score), ]

tail(gene_matrix)
```

<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>activity_score</th><th scope=col>cluster</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1304</th><td>0.164863523405036</td><td>0 </td></tr>
	<tr><th scope=row>510</th><td>0.175571072864751</td><td>4 </td></tr>
	<tr><th scope=row>533</th><td>0.194898488140979</td><td>4 </td></tr>
	<tr><th scope=row>10960</th><td>0.207488050280486</td><td>0 </td></tr>
	<tr><th scope=row>10185</th><td>0.236051070783618</td><td>0 </td></tr>
	<tr><th scope=row>940</th><td>0.317189228896522</td><td>15</td></tr>
</tbody>
</table>



This code chunk displays the cluster that was found to have the highest gene activity score. 


```R
cluster_max <- gene_matrix[gene_matrix$activity_score == max(gene_matrix$activity_score), ]
cluster <- cluster_max$cluster
cluster
```


'15'


This gene activity score can be seen as an extention of research question 1. 
While the data frame ***cluster_gene_dff*** just displays if the clusters differ within their co-accessbile regions, this data frame provides information about the likelihood of a cluster to express a gene. 

Cicero links open promoter regions with the co-accessibility of distal sites. It remains unclear if a data base was used to connect distal sites with known transcription factors. 

This is an approach for further investigation of the research question 1: **Exists an open promoter region in two cell types with different accessibility of distal elements?**. 
With the help of the data frame the differences in the clusters can now be verified or negated. 

The clusters that have gene activitiy scores for the gene of choice can be displayed. If the gene is found in a cluster that has already been displayed (in contrast to the other cluster) by the data frame with the connections, it can be a hint that the peak that differs is decisive when it comes to gene expression. 

***
## Comparison of Promoters in over cell types (Research Question 2) 

Basend on the annotation of peaks to promoter regions, an interesting question would be if a promoter that was previously predicted to be open and active in terms of gene expression is only found in one tissue. 

The second biological question `Cicero` could answer is then:
### Is there a specific promoter region that is only found in one tissue type, comparing all tissues with one another. 


CAUTION: just run this code chunk if the tissue isn't in ***conns_all*** yet. 

Add the column ***tissue***. 
This data frame is now filtered to display the values of ***gene*** just once.
Store the variable ***gene_dff*** as .csv in the directory `/mnt/workspace_stud/stud10/output/`.


```R
# paste the tissue name of the tissue you are working on.
dff$tissue <- "esophagus_mucosa"

#remove duplicates
genes_dff <- dff[!duplicated(dff$gene),]
head(genes_dff)

# append conns_all.csv
write.table(genes_dff, '/mnt/workspace_stud/stud10/output/conns_all.csv', sep = ',', append = TRUE)
```


<table class="dataframe">
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th><th scope=col>tissue</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>138115</th><td>chr17_35766915_35767315 </td><td>chr17_35889180_35889580 </td><td>-0.03748543</td><td>7SK      </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>119972</th><td>chr12_9214539_9214939   </td><td>chr12_9115760_9116160   </td><td>-0.08944022</td><td>A2M      </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>118707</th><td>chr12_8831963_8832363   </td><td>chr12_8830944_8831344   </td><td>-0.02221027</td><td>A2ML1-AS1</td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>438990</th><td>chr3_152243247_152243647</td><td>chr3_151813768_151814168</td><td>-0.03198780</td><td>AADAC    </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>438994</th><td>chr3_152244141_152244541</td><td>chr3_151773081_151773481</td><td>-0.09276838</td><td>AADACP1  </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>170155</th><td>chr17_68782824_68783224 </td><td>chr17_69060849_69061249 </td><td> 0.21821901</td><td>ABCA9    </td><td>esophagus_mucosa</td></tr>
</tbody>
</table>



Input the data frame which stores all active promoters per tissue. 


```R
conns_all <- read.csv("/mnt/workspace_stud/stud10/output/conns_all.csv", 
                      sep = ',', row.names = NULL)
conns_all <- conns_all[c("Peak1", "Peak2", "coaccess", "gene", "tissue")]
head(conns_all)
```


<table class="dataframe">
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th><th scope=col>tissue</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr17_35766915_35767315 </td><td>chr17_35889180_35889580 </td><td>-0.0374854318970528</td><td>7SK      </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>2</th><td>chr12_9214539_9214939   </td><td>chr12_9115760_9116160   </td><td>-0.089440218142176 </td><td>A2M      </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>3</th><td>chr12_8831963_8832363   </td><td>chr12_8830944_8831344   </td><td>-0.0222102706602923</td><td>A2ML1-AS1</td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>4</th><td>chr3_152243247_152243647</td><td>chr3_151813768_151814168</td><td>-0.031987797677086 </td><td>AADAC    </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>5</th><td>chr3_152244141_152244541</td><td>chr3_151773081_151773481</td><td>-0.0927683812435143</td><td>AADACP1  </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>6</th><td>chr17_68782824_68783224 </td><td>chr17_69060849_69061249 </td><td>0.218219011600414  </td><td>ABCA9    </td><td>esophagus_mucosa</td></tr>
</tbody>
</table>



Check if there is a unique value in the ***gene*** column. In order to do so the data frame has to be filtered by the counts of the column ***gene***, this is done by `counts()`. The filter then is only to take genes that the count-table assigns with 1. In the end the variable ***f*** stores a data frame with one column ***gene***. 


```R
f <- conns_all %>% 
    count(gene) %>% 
    filter(n == 1) %>% 
    select(-n)
```

The data frame with all genes is now modified to display only rows which value in the column ***gene*** is found in the data frame ***f*** by `inner_join()`. 


```R
unique_genes <- inner_join(conns_all, f, by = "gene")
unique_genes <- unique_genes[, c("Peak2", "gene", "tissue")]

unique_genes[order(unique_genes$gene), ]
```


<table class="dataframe">
<caption>A data.frame: 2783 × 3</caption>
<thead>
	<tr><th></th><th scope=col>Peak2</th><th scope=col>gene</th><th scope=col>tissue</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1808</th><td>chr19_58347525_58347925  </td><td>A1BG-AS1 </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>618</th><td>chr12_8844678_8845078    </td><td>A2ML1    </td><td>lung            </td></tr>
	<tr><th scope=row>1</th><td>chr12_8830944_8831344    </td><td>A2ML1-AS1</td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>1346</th><td>chr12_125092361_125092761</td><td>AACS     </td><td>colon_transverse</td></tr>
	<tr><th scope=row>2</th><td>chr3_151813768_151814168 </td><td>AADAC    </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>1809</th><td>chr17_76453024_76453424  </td><td>AANAT    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1810</th><td>chr6_44313155_44313555   </td><td>AARS2    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>619</th><td>chr17_69175705_69176105  </td><td>ABCA10   </td><td>lung            </td></tr>
	<tr><th scope=row>620</th><td>chr1_94021241_94021641   </td><td>ABCA4    </td><td>lung            </td></tr>
	<tr><th scope=row>1347</th><td>chr19_1045204_1045604    </td><td>ABCA7    </td><td>colon_transverse</td></tr>
	<tr><th scope=row>3</th><td>chr2_168926976_168927376 </td><td>ABCB11   </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>1348</th><td>chr12_122944570_122944970</td><td>ABCB9    </td><td>colon_transverse</td></tr>
	<tr><th scope=row>1811</th><td>chr17_50634626_50635026  </td><td>ABCC3    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1812</th><td>chr1_94418634_94419034   </td><td>ABCD3    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>621</th><td>chr7_151226687_151227087 </td><td>ABCF2    </td><td>lung            </td></tr>
	<tr><th scope=row>4</th><td>chr2_43831669_43832069   </td><td>ABCG8    </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>1813</th><td>chr2_27123563_27123963   </td><td>ABHD1    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1814</th><td>chr3_100993271_100993671 </td><td>ABI3BP   </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1815</th><td>chr6_139028452_139028852 </td><td>ABRACL   </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>5</th><td>chr17_37088906_37089306  </td><td>ACACA    </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>622</th><td>chr2_210224465_210224865 </td><td>ACADL    </td><td>lung            </td></tr>
	<tr><th scope=row>1816</th><td>chr3_195442199_195442599 </td><td>ACAP2    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>623</th><td>chr11_108122155_108122555</td><td>ACAT1    </td><td>lung            </td></tr>
	<tr><th scope=row>624</th><td>chr17_63477205_63477605  </td><td>ACE      </td><td>lung            </td></tr>
	<tr><th scope=row>1817</th><td>chr17_75978714_75979114  </td><td>ACOX1    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>625</th><td>chr2_110732762_110733162 </td><td>ACOXL    </td><td>lung            </td></tr>
	<tr><th scope=row>626</th><td>chr19_11577577_11577977  </td><td>ACP5     </td><td>lung            </td></tr>
	<tr><th scope=row>6</th><td>chr19_6161385_6161785    </td><td>ACSBG2   </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>627</th><td>chr7_5527991_5528391     </td><td>ACTB     </td><td>lung            </td></tr>
	<tr><th scope=row>1818</th><td>chr17_81512262_81512662  </td><td>ACTG1    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>2768</th><td>chr1_247007912_247008312</td><td>ZNF695     </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2769</th><td>chr19_52993358_52993758 </td><td>ZNF702P    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2770</th><td>chr19_56595047_56595447 </td><td>ZNF71      </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2771</th><td>chr3_75784983_75785383  </td><td>ZNF717     </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2772</th><td>chr19_23914637_23915037 </td><td>ZNF726     </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1804</th><td>chr19_23003022_23003422 </td><td>ZNF728     </td><td>colon_transverse</td></tr>
	<tr><th scope=row>611</th><td>chr7_64307879_64308279  </td><td>ZNF736     </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>2773</th><td>chr6_35258801_35259201  </td><td>ZNF76      </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1805</th><td>chr7_150382500_150382900</td><td>ZNF775     </td><td>colon_transverse</td></tr>
	<tr><th scope=row>612</th><td>chr7_150383848_150384248</td><td>ZNF775-AS1 </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>2774</th><td>chr19_12092020_12092420 </td><td>ZNF788P    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1343</th><td>chr19_34964084_34964484 </td><td>ZNF792     </td><td>lung            </td></tr>
	<tr><th scope=row>1344</th><td>chr19_12400634_12401034 </td><td>ZNF799     </td><td>lung            </td></tr>
	<tr><th scope=row>2775</th><td>chr2_184598266_184598666</td><td>ZNF804A    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1806</th><td>chr7_88759317_88759717  </td><td>ZNF804B    </td><td>colon_transverse</td></tr>
	<tr><th scope=row>2776</th><td>chr19_20424788_20425188 </td><td>ZNF826P    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2777</th><td>chr19_52171480_52171880 </td><td>ZNF836     </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2778</th><td>chr19_9792953_9793353   </td><td>ZNF846     </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>1807</th><td>chr19_22784000_22784400 </td><td>ZNF99      </td><td>colon_transverse</td></tr>
	<tr><th scope=row>2779</th><td>chr6_30060612_30061012  </td><td>ZNRD1ASP   </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>613</th><td>chr7_50093073_50093473  </td><td>ZPBP       </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>614</th><td>chr17_39867890_39868290 </td><td>ZPBP2      </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>615</th><td>chr1_26225232_26225632  </td><td>ZPLD2P     </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>616</th><td>chr19_58033731_58034131 </td><td>ZSCAN1     </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>2780</th><td>chr6_28399605_28400005  </td><td>ZSCAN12    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2781</th><td>chr6_28124365_28124765  </td><td>ZSCAN16    </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2782</th><td>chr6_28137151_28137551  </td><td>ZSCAN16-AS1</td><td>artery_tibial   </td></tr>
	<tr><th scope=row>2783</th><td>chr6_116668613_116669013</td><td>ZUP1       </td><td>artery_tibial   </td></tr>
	<tr><th scope=row>617</th><td>chr17_4074831_4075231   </td><td>ZZEF1      </td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>1345</th><td>chr1_77631748_77632148  </td><td>ZZZ3       </td><td>lung            </td></tr>
</tbody>
</table>




```R
len <- unique_genes[unique_genes$tissue == 'esophagus_mucosa', ]
length(len$tissue)
```


617


This data frame is providing information about uniquely found promoter regions over all considered tissues. The column ***Peak2*** is the genomic location of the promoter region of the gene, whose name is displayed in ***gene***, ***tissue*** tells in which tissue this unique promoter region was found. 

The initial question was **is there a specific promoter region that is only found in one tissue type, comparing all tissues with one another.**

The data frame displays: there are indeed promoters to genes that are uniquely found in one tissue. A way to further investigate here is to have a closer look on the function the gene holds. This could answer the question why it is only accessible via one single tissue. 
One has to be careful, because Cicero only gives hints and predictions. What should follow are laboratory studies that further support the thesis of gene expression, this could be done by PCR or an RNA-seq assay.

***
## Cicero Connection Plot


```R
library(cicero)
library(ggplot2)
library(stringr)
```

The Cicero package includes a plotting functing. The function `plot_connections()` visualizes the pairwise co-accessibility. By usage of the `Gviz` framework `Cicero`plots genome browser-style plots. This function has many options. The function first requires the connection-score data frame, a frame of connections including the columns ***Peak1***, ***Peak2*** and ***coaccess***. This information is provided by the data frame ***conns***. Furthermore a chromosome on which the plotting has to be performed must be given as well as the information on ***minbp*** and ***maxbp***, base pair coordinates of start and end region of the plot.
To automate the plotting variables are used, which are likely to show good plot results. 

First the GTF File is loaded as ***gene_anno*** variable, some columns are renamed to match requirements. 


```R
gene_anno <- rtracklayer::readGFF("/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.gtf")
# rename some columns to match requirements
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name
```

In order to create the variables ***minbp***, ***maxbp*** and ***chrom*** a subset of the ***dff*** data frame is used. The ***cutoff*** of `plot_connections()` is set to 0.2, lower scores are not plottet, so this will be the threshold here. 


```R
filtered_dff <- subset(dff, coaccess > 0.2)
```

A gene can be chosen or the variable ***g*** can be used, which displays the ***gene*** which is representing the gene that has the most co-accessbility scores. The last 10 entries will be taken into consideration. 


```R
g <- sort(table(filtered_dff$gene))
g <- tail(g, n=1)
g <- names(g)
dff_gene <- filtered_dff[filtered_dff$gene == g,]
dff_gene <- tail(dff_gene, n = 10)
tail(dff_gene)
```


<table class="dataframe">
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th><th scope=col>gene</th><th scope=col>tissue</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>385707</th><td>chr2_219408563_219408963</td><td>chr2_219451703_219452103</td><td>0.6559549</td><td>SPEG</td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>385959</th><td>chr2_219453711_219454111</td><td>chr2_219460586_219460986</td><td>0.6661128</td><td>SPEG</td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>385708</th><td>chr2_219408563_219408963</td><td>chr2_219460586_219460986</td><td>0.6991157</td><td>SPEG</td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>386007</th><td>chr2_219460586_219460986</td><td>chr2_219461011_219461411</td><td>0.7084824</td><td>SPEG</td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>386022</th><td>chr2_219461011_219461411</td><td>chr2_219460586_219460986</td><td>0.7084824</td><td>SPEG</td><td>esophagus_mucosa</td></tr>
	<tr><th scope=row>385802</th><td>chr2_219433809_219434209</td><td>chr2_219460586_219460986</td><td>0.7612388</td><td>SPEG</td><td>esophagus_mucosa</td></tr>
</tbody>
</table>



To now extract the base pair coordinates the data frame has to be rearranged and the columns has to be splittet. After `order()` the smallest ***Peak1*** is displayed and can be stored in ***minbp*** by splitting the column ***Peak1*** and taking the first entry. ***chrom*** is initialized with the newly created column ***chr***, while  ***view*** is initialized with the ***Peak2*** column. 


```R
minbp <- dff_gene[order(dff_gene$Peak1), ]
minbp <- minbp[1, ]
minbp[c('chr', 'bp1', 'bp2')] <- str_split_fixed(minbp$Peak1, '_', 3)
view <- minbp$Peak2
chrom <- minbp$chr
minbp <- as.numeric(minbp$bp1)
```

This `order()` process is to be repeated for ***maxbp*** only taking ***Peak2*** into account and choosing the last entry.


```R
maxbp <- dff_gene[order(dff_gene$Peak2), ]
maxbp <- tail(maxbp, n = 1 )
maxbp[c('chr', 'bp1', 'bp2')] <- str_split_fixed(maxbp$Peak2, '_', 3)
maxbp <- as.numeric(maxbp$bp2)
```

The function `plot_connections()` has the options, besides the input data frame and the coordinates to plot on, to set a ***viewpoint***, which sets a line to see connections originating from a specific place in the genome, e.g. the base pair coordinates of a promoter region. The ***gene_model*** is the reference GTF data frame, in this project ***gene_anno***. 

Further Options that seem to be important for this workpackage are: 
The ***coaccess_cutoff*** to define the minimum co-accessibility score that should be plottet, the ***connection_width*** to define the width of the connection lines. Moreover the `Cicero` website recommends to set ***collapseTranscripts*** to ***longest*** to determine how and whether to collaps related transcripts. 


```R
plot_connections(dff, chrom, minbp, maxbp,
                 viewpoint = view,
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.50, 
                 connection_width = 1.5, 
                 collapseTranscripts = "longest" )
```


![output_75_0](https://user-images.githubusercontent.com/93346891/160593455-b04a4737-aa63-4cab-91a0-028b2c270116.png)


Visualized by this plot are the connections from a peak that is assigned to a promoter region in the genome.
The red bars represents the peaks, under there the genome window is displayed as grey arrows. Is a gene assigned to this genome window it is shown in the green bar with the corresponding name in front. 

Many distal peaks are found, some on other genes and some in the genome distal from the assigned promoter region. 
