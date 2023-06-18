set -o pipefail
set -o errexit

declare MD_BARCODE_TAG
if [ ! -z "$6" ]; then
  MD_BARCODE_TAG="BARCODE_TAG=$6"
/usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=$5 METRICS_FILE=$4 QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT "$MD_BARCODE_TAG" | /usr/bin/sambamba sort -t $2 -m 18G -o $3 /dev/stdin
else
  /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=$5 METRICS_FILE=$4 QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | /usr/bin/sambamba sort -t $2 -m 18G -o $3 /dev/stdin
fi
