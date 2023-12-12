##cat GSE87371_series_matrix.txt |grep -v "\!" > GSE87371_series_matrix_no_anno.txt
##cat GPL570-55999.txt|grep -v "#" > GPL570-55999_no_anno.txt

##读取数据信息
setwd("E://曦昀项目//2.双疾病//1.BRCA//1.数据下载")
expr_df <- read.table("GSE11783.txt",header=TRUE) 
class(expr_df) # dataframe
dim(expr_df)
expr_df[1:3,]

sample_names <- colnames(expr_df)[-1]
probe_ids <- expr_df$ID_REF

##读取平台信息
id_table <- read.table("GPL570-55999.txt",header=TRUE,sep = "\t",dec = ".",comment.char = "#",fill = T,quote = "")
id_table[1:11,1:11]
colnames(id_table)

#install.packages("data.table")
require(data.table)

probe2symbol <- id_table[,c("ID","Gene.Symbol")]

# 得到基因symbol，一个探针可能对应多个基因。
symbol <-tstrsplit(id_table$Gene.Symbol, "///", fixed=TRUE)[[1]]
# 去除空格
symbol<- trimws(symbol, which = c("both", "left", "right"),whitespace = "[ \t\r\n]")
probe2symbol["symbol"] <- symbol
# 去掉gene_assignment列
probe2symbol <-  probe2symbol[,c("ID","symbol")]

head(probe2symbol)
# ID列名改为probe_id
colnames(probe2symbol) <- c("ID_REF","symbol")



#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("hugene10sttranscriptcluster.db")

library(hugene10sttranscriptcluster.db)

## Bimap interface:
x <- hugene10sttranscriptclusterSYMBOL

# Get the probe identifiers that are mapped to a gene symbol 
mapped_probes <- mappedkeys(x)
# Convert to dataframe
probe2symbol2 <- as.data.frame(x[mapped_probes])

# probe2symbol2[which(probe2symbol2$probe_id == "7898916"),]

# check
dim(probe2symbol2)
probe2symbol[1:30,]

probe2symbol2[which(probe2symbol2$probe_id == "7896859"),]
probe2symbol[which(probe2symbol2$probe_id == "7896859"),]
id_table[which(probe2symbol2$probe_id == "7896859"),c("Gene.Symbol")]


merged_expr_df <- merge(x =probe2symbol , y = expr_df)

# 去掉probe_id没有对应基因symbol的行
filt_expr_df <- merged_expr_df[complete.cases(merged_expr_df),]
#针对某一列过滤，本例效果一样。
filt_expr_df <- merged_expr_df[complete.cases(merged_expr_df[,c("symbol")]),]
# check
dim(merged_expr_df)
dim(filt_expr_df)

table(filt_expr_df$symbol) # 有重复，多个探针对应一个基因，不能作为行名

filt_expr_df[1:3,]

# 去掉ID_REF列(probe id)
filt_expr_df <- subset(filt_expr_df, select = -ID_REF)
# 取每个基因所有探针的平均值或最大值作为基因的表达量
m_df <- aggregate(.~symbol,data=filt_expr_df,mean)
m_df <- aggregate(.~symbol,data=filt_expr_df,max)
dim(m_df)
# 查看结果
filt_expr_df[which(filt_expr_df$symbol == "ARHGDIA"),]
m_df[which(m_df$symbol == "ARHGDIA"),]

rownames(m_df) <- m_df$symbol

m_df <- subset(m_df, select = -symbol) #去掉symbol列
exprSet <- as.matrix(m_df)
head(exprSet)

write.table(exprSet,"GEOmatrix.txt",sep="\t",quote = FALSE)
