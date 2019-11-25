from ftplib import FTP
import subprocess

site = "ftp.1000genomes.ebi.ac.uk"
ftp_directory = "vol1/ftp/phase3/data/HG00096/alignment/"
bam_file = "HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
# gatk_path = "../../../Downloads/gatk-4.1.4.0/gatk-4.1.2.0/gatk"

def ftp_file_download(ftp_site, directory, file_list):
    # download to current directory; can add support for checking if sufficient space for download to current directory
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd(ftp_directory)
    # ftp.retrlines("LIST")
    for filename in file_list:
        file_ptr = open(filename, "wb")
        ftp.retrbinary("RETR " + filename, file_ptr.write)
        file_ptr.close()


def subprocess_exec(cmd):
    print("Run command: {}".format(cmd))
    p1 = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p1.communicate()
    if out.decode() != '':
        print("\nOUTPUT\n{}".format(out.decode()))
    if err.decode() != '':
        print("\nMESSAGE\n{}".format(err.decode()))


def main():

    # download files
    files_to_download = [bam_file, bam_file + '.bai']  # bam and bai files
    ftp_file_download(site, ftp_directory, files_to_download)

    # did a quick view of header samtools view -H and first few lines on terminal samtools view <filename> | less

    # sort - in this case it's already sorted
    cmd = """samtools sort {} > {}.sorted.bam""".format(files_to_download[0], files_to_download[0][0:-4])
    subprocess_exec(cmd)

    # observe bam stats
    cmd = """samtools flagstat {}.sorted.bam""".format(files_to_download[0][0:-4])
    subprocess_exec(cmd)

    # get insert size metrics with Picard CollectInsertMetrics
    # ASSUMPTION 1: fragment size approx equals insert size
    # here i ignore adapter length since actually fragment size = insert size + length of adapters at both ends
    # ASSUMPTION 2: insert size approx equals TLEN ( field 9 of sam ) which I believe Picard CollectInsertMetrics uses
    # this assumption would be violated in case of RNAseq data due to splicing and hence introns not being sequenced
    cmd = """gatk CollectInsertSizeMetrics -INPUT {}.sorted.bam -Histogram_FILE insert_size_distribution_plot.pdf \
                -OUTPUT insert_size_information.txt""".format(files_to_download[0][0:-4])
    subprocess_exec(cmd)


if __name__ == "__main__":
    main()


