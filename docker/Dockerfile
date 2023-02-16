from docker.io/python:3.6
RUN pip install pyvcf -i https://mirrors.aliyun.com/pypi/simple/ --trusted-host mirrors.aliyun.com
COPY openanno /usr/local/bin/openanno
COPY vcf2annoinput_somatic.py /usr/local/bin/vcf2annoinput_somatic.py
COPY vcf2annoinput_germline.py /usr/local/bin/vcf2annoinput_germline.py
RUN chmod a+x /usr/local/bin/openanno /usr/local/bin/vcf2annoinput_somatic.py /usr/local/bin/vcf2annoinput_germline.py
