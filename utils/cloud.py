from azure.common.client_factory import get_client_from_cli_profile
from azure.mgmt.compute import ComputeManagementClient
from azure.mgmt.network import NetworkManagementClient
import os, uuid, sys
from azure.storage.blob import BlockBlobService, PublicAccess
import subprocess
import json
import tarfile
import shutil
import gzip

class ResultsInterface(object):

    def __init__(self, version="", dir="./tmp"):
        self.dir = dir + version
        try:
            os.makedirs(dir)
        except Exception as e:
            pass
        self.results = []
        self.version = version
        cmd = "az storage blob list --container-name cellranger{} --account-name scrnadata --auth-mode login".format(version)
        json_dump = subprocess.check_output(cmd.split())
        blobs = json.loads(json_dump)
        for blob in blobs:
            if "patient1" not in blob["name"]: continue
            result = TenxDataStorage(blob["name"], version=version)
            result.download()
            self.results.append(result)


class TenxDataStorage(object):

    def __init__(self, sampleid, version="v3"):
        self.sampleid = sampleid.replace(".tar.gz","")
        self.storage_account = "scrnadata"
        self.container = "cellranger{}".format(version)
        self.rdatacontainer = "rawrdata{}".format(version)
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token='?sv=2018-03-28&ss=bfqt&srt=sco&sp=rwdlacup&se=2021-03-19T02:52:48Z&st=2019-02-22T19:52:48Z&spr=https&sig=4oAGvIyqi9LPY89t21iWWp4XbbIgpmaldgcwWUOuf14%3D')
        self.tenx_path = None
        self.cache = ".cache"
        try:
            os.makedirs(self.cache)
        except Exception as e:
            pass

    def rdata(self):
        local = ".cache/{}.rdata".format(self.sampleid)
        raw = "{}.rdata".format(self.sampleid)
        if not os.path.exists(local):
            self.block_blob_service.get_blob_to_path(self.rdatacontainer, raw, local)
        return local

    def upload_cellranger(self, tenx):
        bam_tarball = tenx.bam_tarball()
        bam_tarball_name = os.path.split(bam_tarball)[-1]
        outs_tarball = tenx.outs_tarball()
        outs_tarball_name = os.path.split(outs_tarball)[-1]
        bam_blob_path = "{sampleid}/{local}".format(self.sampleid, bam_tarball_name)
        self.upload("bams" ,bam_blob_path, bam_tarball)
        outs_blob_path = "{sampleid}/{local}".format(self.sampleid, bam_tarball_name)
        self.upload("cellranger", outs_blob_path, bam_tarball)


    def unpack(self,path):
        tar = tarfile.open(path)
        tar.extractall(".cache")
        tar.close()
        assert os.path.exists(".cache/outs"), "Failed to download."
        os.rename(".cache/outs",os.path.join(self.cache,self.sampleid))


    def download(self):
        local = ".cache/{}.tar.gz".format(self.sampleid)
        gzipped = "{}.tar.gz".format(self.sampleid)
        if not os.path.exists(local):
            self.block_blob_service.get_blob_to_path(self.container, gzipped, local)
            self.unpack(local)
        self.tenx_path = os.path.join(self.cache,self.sampleid)
        return self.tenx_path

    def upload(self,container, blob, local):
        print ("Uploading {} to {} in {}".format(container,local,blob))
        cmd = "az storage blob upload --container-name {} --file {} --name {} --auth-mode login".format(container,local,blob)
        subprocess.call(cmd.split())

class FastqDataStorage(object):

    def __init__(self, sampleid):
        self.sampleid = sampleid
        self.storage_account = "scrnadata"
        self.container = "rnaseq"
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token='?sv=2018-03-28&ss=bfqt&srt=sco&sp=rwdlacup&se=2021-03-19T02:52:48Z&st=2019-02-22T19:52:48Z&spr=https&sig=4oAGvIyqi9LPY89t21iWWp4XbbIgpmaldgcwWUOuf14%3D')

    def set_data_path(self, path):
        self.datapath = path
        try:
            os.makedirs(path)
        except Exception as e:
            pass

    def download_fastqs(self):
        fastqs = self.block_blob_service.list_blobs(self.container)
        for fastq in fastqs:
            if len(fastq.name.split("/")) != 2: continue
            sample, fastqname = fastq.name.split("/")
            if sample == self.sampleid:
                local = os.path.join(self.datapath, fastqname)
                if not os.path.exists(local):
                    self.block_blob_service.get_blob_to_path(self.container, fastq.name, local)
        return self.datapath

class ReferenceDataStorage(object):

    def __init__(self,build, referencepath):
        self.sampleid = build
        self.storage_account = "scrnadata"
        self.container = "reference"
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token='?sv=2018-03-28&ss=bfqt&srt=sco&sp=rwdlacup&se=2021-03-19T02:52:48Z&st=2019-02-22T19:52:48Z&spr=https&sig=4oAGvIyqi9LPY89t21iWWp4XbbIgpmaldgcwWUOuf14%3D')
        self.reference = os.path.join(referencepath, build)
        self.build = build
        self.referencepath = referencepath

    def extract(self,path):
        tar = tarfile.open(path)
        tar.extractall(path=self.referencepath)
        tar.close()


    def download(self):
        if not os.path.exists(self.reference):
            local = os.path.join(self.referencepath, "{}.tar.gz".format(self.build))
            self.block_blob_service.get_blob_to_path(self.container, "{}.tar.gz".format(self.build), local)
            self.extract(local)
        return self.reference


class VirtualMachine(object):

    def __init__(self):
        compute_client = get_client_from_cli_profile(ComputeManagementClient)
        network_client = get_client_from_cli_profile(NetworkManagementClient)
        nic = network_client.network_interfaces.get(GROUP_NAME, NIC_NAME)
        print('\nCreating Linux Virtual Machine')
        vm_parameters = create_vm_parameters(nic.id, VM_REFERENCE['linux'])
        print(vm_parameters)
        async_vm_creation = compute_client.virtual_machines.create_or_update(GROUP_NAME, VM_NAME, vm_parameters)
        async_vm_creation.wait()

    def run_job(self):
        pass

    def check_status(self):
        pass

if __name__ == '__main__':
    test_fastqs = FastqDataStorage("test-fastqs")
    test_fastqs.download_fastqs()
