#! /bin/bash

go build -o docker/openanno main.go

cd $(dirname $0)

version=$(head -1 version)
sub_version1=$(echo $version|awk -F '.' '{print $NF}')
sub_version2=$(($sub_version1+1))
version=$(echo $version|sed 's:'$sub_version1'$:'$sub_version2':g')

tag=registry.cn-shanghai.aliyuncs.com/kszy-biosoft/openanno:$version

docker build -t $tag . && docker push $tag && echo $version >version 
