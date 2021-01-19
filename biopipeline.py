from pathlib import Path
import datetime
import subprocess


class jmpipe:
    """Controler of pipeline

    """

    def __init__(self, manifest, name='RNA-seq', rawdata='rawdata', cleandata: str = './cleandata',
                 refdir: str = './ref', alined: str = './alined'):
        # software
        self.fastp = '/home/mwshi/jimmy/software/fastp'
        self.multiqc_ = '/home/mwshi/miniconda2/bin/multiqc'
        self.fastqc = '/home/mwshi/jimmy/software/FastQC/fastqc'
        self.hisat2 = '/home/zxchen/anaconda3/bin/hisat2'
        self.samtools = '/home/mwshi/miniconda2/bin/samtools'

        # manifest
        self.manifest = manifest
        self.sampledict = self.getrawdata()
        self._name = name

        # info
        self._time = datetime.datetime.now()

        # path
        self.mypwd = Path.cwd()
        self.cleandata = cleandata  # after QC
        self.refdir = refdir  # reffile
        self.alined = alined  # afteralined:bam|bam.fai
        self.rawdata = rawdata
        self.makedir = self.makedir()

    def makedir(self):
        for i in [self.refdir, self.cleandata, self.alined, self.rawdata]:
            p = Path(i)
            if not p.exists():
                p.mkdir()
                print('===>made dir\t:{}'.format(i))
        return

    def getrawdata(self):
        sampdict = {}
        with open(self.manifest, 'r') as f:
            for i in f:
                datadict = {}
                sample_id, datadir, f, r = i.strip().split()
                f, r = Path(datadir) / Path(f), Path(datadir) / Path(r)
                datadict['rawfq'] = (str(f), str(r))
                sampdict[sample_id] = datadict
        return sampdict

    @property
    def samples(self):
        samlelist = []
        for sample_id, v in self.sampledict.items():
            f, r = list(v.values())[0]
            samlelist.append('sample_id:{}:\tForword:{}\tReverse:{}'.format(sample_id, f, r))
        return samlelist

    @property
    def info(self):
        return 'The project is ===>{}\tmade in {} by jmzhang'.format(self._name, self._time)


class RnaSeq(jmpipe):
    """This is a pipeline of RNA-Seq
    step01: project = RnaSeq('pathtomanifest.txt', 'projectname')
            project.info -> project-name project-time
            project.samples -> manifest.txt
    step02: QC:
                => project.qc(rawqc=True, cleanqc=True, factp=True)
                => project.mutiqc()
    """

    def qc(self, rawqc=True, cleanqc=True, factp=True):
        def fastqc(filelist: list, outpath):
            for _ in filelist:
                cmd = '{} -o {}/fastqc_res -f fastq {}'.format(self.fastqc, outpath, _)
                print(cmd)

        if rawqc:
            """fastqc for raw data
            """
            p = Path(self.rawdata) / 'rawqc'
            if not p.exists():
                p.mkdir()
            fastqlist = []
            for sample_id, _ in self.sampledict.items():
                f, r = _['rawfq']
                fastqlist.extend([f, r])
            fastqc(fastqlist, p)

        if factp:
            for sample_id, _ in self.sampledict.items():
                f, r = _['rawfq']
                fout, rout = f.replace('fastq', 'fastp'), r.replace('fastq', 'fastp')
                _['cleanfq'] = (fout, rout)
                cmd = '{} -i {} -I {} -o {} -O {}'.format(self.fastp, f, r, fout, rout)
                print(cmd)

        if cleanqc:
            leanfastq = []
            for x in self.sampledict.values():
                f, r = x['cleanfq']
                leanfastq.extend([f, r])
            fastqc(leanfastq, self.cleandata)

        pass

    def multiqc(self,cleanfq=True,rawfq=True):
        if cleanfq:
            cmd = '{} {} -O {}'.format(self.multiqc_,str(self.mypwd / 'cleandata/fastqc_res'), str(self.mypwd / 'multiqc_res'))
            print(cmd)
        if rawfq:
            cmd = '{} {} -O {}'.format(self.multiqc_,str(self.mypwd / 'rawdata/fastqc_res'), str(self.mypwd / 'multiqc_res'))
            print(cmd)

    def align(self):

        pass

    def getre(self):

        pass


manifest = r'E:\repository\pycharm_project\manifest'

project1 = RnaSeq(manifest, 'crc-lncRNA')
# print(project1.samples)
# print(project1.info)

project1.qc(rawqc=True, factp=True, cleanqc=True)
print('='*30)
#print(project1.sampledict)
project1.multiqc()