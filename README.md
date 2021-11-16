# OpenAnno

#### 介绍

OpenAnno
是一款基于GO语言开发的Gene注释软件，它受[ANNOVAR](https://annovar.openbioinformatics.org/en/latest/)启发，根据作者时间工作经验需求因现存VEP、ANNOVAR等注释软件不能完全满足业务需求而重写，力求做到代码架构简洁，可修改性高。

#### 安装教程

```shell
git clone https://gitee.com/zhuying/open-anno
cd open-anno
go build -o openanno main.go
```

#### 使用说明

##### 数据准备

1. refGene.txt
2. ref.fa
3. humandb: directory of database
4. database annotation file: such as gnomad.txt
5.

```shell
samtools dict ref.fa > humandb/ref.dict
openanno prepare transcript -d humandb/ -r ref.fa -g refGene.txt -o humandb/refgene
openanno prepare database -i gnomad.txt -o humandb/gnomad
```

6. 修改humandb/config.yaml

##### 运行注释

```shell
openanno anno snv -d humandb/ -i test.avinput -o test.json
```

#### 参与贡献

1. Fork 本仓库
2. 新建 Feat_xxx 分支
3. 提交代码
4. 新建 Pull Request