from pathlib import Path
import datetime
import subprocess


class JmPipe:
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


class RnaSeq(JmPipe):
    """This is a pipeline of RNA-Seq
    step01: project = RnaSeq('pathtomanifest.txt', 'projectname')
            project.info -> project-name project-time
            project.samples -> manifest.txt
    step02: QC:
                => project.qc(rawqc=True, cleanqc=True, factp=True)
                => project.mutiqc()
    """

    def qc(self, rawqc=True, cleanqc=True, fastp=True):
        def fastqc(filelist: list, outpath:Path):
            for _ in filelist:
                if not (outpath / 'fastqc_res').exists():
                    (outpath / 'fastqc_res').mkdir()
                cmd = '{} -o {}/fastqc_res -f fastq {}'.format(self.fastqc, outpath, _)
                subprocess.run(cmd, shell=True)

        if rawqc:
            """fastqc for raw data
            """
            fastqlist = []
            for sample_id, _ in self.sampledict.items():
                f, r = _['rawfq']
                fastqlist.extend([f, r])
            fastqc(fastqlist, Path(self.rawdata))

        if fastp:
            for sample_id, _ in self.sampledict.items():
                f, r = _['rawfq']
                fn, rn = Path(f).name, Path(r).name
                fout, rout = fn.replace('fastq', 'fastp'), rn.replace('fastq', 'fastp')
                fout, rout = Path(self.mypwd)/ 'cleandata' / fout, Path(self.mypwd) / 'cleandata' / rout
                _['cleanfq'] = (fout, rout)
                cmd = '{} -i {} -I {} -o {} -O {}'.format(self.fastp, f, r, fout, rout)
                print(cmd)
                #subprocess.run(cmd, shell=True)

        if cleanqc:
            leanfastq = []
            for x in self.sampledict.values():
                f, r = x['cleanfq']
                leanfastq.extend([f, r])
            fastqc(leanfastq, Path(self.cleandata))

    def multiqc(self,cleanfq=True,rawfq=True):
        if cleanfq:
            cmd = '{} {} -O {}'.format(self.multiqc_,str(self.mypwd / 'cleandata/fastqc_res'), str(self.mypwd / 'multiqc_res'))
            subprocess.run(cmd, shell=True)
        if rawfq:
            cmd = '{} {} -O {}'.format(self.multiqc_,str(self.mypwd / 'rawdata/fastqc_res'), str(self.mypwd / 'multiqc_res'))
            subprocess.run(cmd, shell=True)

    def align(self):

        pass

    def getre(self):

        pass


if __name__ ==  '__main__':
    manifest = r'E:\repository\pycharm_project\manifest'

    project1 = RnaSeq(manifest, 'crc-lncRNA')
    # print(project1.samples)
    # print(project1.info)

    project1.qc(rawqc=False, fastp=True, cleanqc=False)
    print('='*30)
    #print(project1.sampledict)
    #project1.multiqc()
    print(project1.sampledict)
