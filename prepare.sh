if [ $3 -ne 3 ]
then 
    echo "bash $0 <input> <output>"
    exit
fi

format=$1
input=$2
output=$3

case format in
gpe)
    zcat $input | sort -k 3,3 -k 5,5n -k 6,6n |bgzip -c > $output
    tabix -s 3 -b 5 -e 6 $output
    ;;
vcf)
    bcftools sort $input -Oz -o $output
    bcftools index --tbi $output
    ;;
bed)
    sort -k 1,1 -k 2,2n -k 3,3n $input |bgzip -c > $output
    tabix -s 1 -b 2 -e 3 $output
    ;;
*)
    echo "error format: $format"
    ;;
esac