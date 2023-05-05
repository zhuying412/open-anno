FROM registry.cn-shanghai.aliyuncs.com/kszy-biosoft/autopvs1:2023-03-31-16-24-48
COPY bin/openanno /usr/local/bin/openanno
RUN chmod a+x /usr/local/bin/openanno
