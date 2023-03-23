from docker.io/tensorflow/tensorflow:2.11.0
RUN pip install pyvcf spliceai -i https://pypi.douban.com/simple/ --trusted-host pypi.douban.com
COPY openanno /usr/local/bin/openanno
COPY spliceai.py /usr/local/bin/spliceai.py
RUN chmod a+x /usr/local/bin/openanno /usr/local/bin/spliceai.py
