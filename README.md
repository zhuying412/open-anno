# OpenAnno

## 介绍

OpenAnno是使用GO语言编写的，用于基因测序分析中SNV、CNV变异注释的软件。本软件以ANNOVAR、VEP等注释软件为灵感，在实现主要注释功能的同时，对生产实践遇到部分问题进行优化调整。

>声明： 本软件以ANNOVAR、VEP注释分析结果为灵感，所有代码均从头编写为参考任何现有代码或算法。

## 拟解决的问题

- SNV注释
  - 基因功能注释
  - 位点数据库注释
  - 区域数据库注释
- CNV注释
  - 基因功能注释
  - 区域数据库注释

## 主要功能

- pre: 数据准备
  - gb: 基因功能注释相关数据准备
  - db: 位点、区域相关数据准备
- anno: 变异注释
  - gb: 基因功能注释
    - snv: SNV基因功能注释
    - cnv: CNV基因功能注释
  - fb: 位点数据库注释，主要针对SNV
  - rb: 区域数据库注释
- tools: 相关工具
  - merge: 合并分析结果，即合并基因功能、位点、区域等多个注释结果到唯一解雇文件
  - vcf2av: VCF转本软件专用数据格式（同ANNOVAR avinput） 
  - av2vcf: 本软件专用数据格式转VCF数据格式

## 优化项

- 现有注释软件核酸与氨基酸改变的[3’端原则问题](https://www.yuque.com/zhuy/bio/aog31h)
- snv deletion 跨过CDS与intron连接区（同时包含CDS与intron部分位点）产生的核酸与氨基酸改变问题