import argparse
import os
from toil.common import Toil
from toil.job import Job

from RepeatAnnotator.repeatAnnotator import addRepeatAnnotatorOptions, readGFF, writeGFF, buildSubfamilies, buildSubfamilies, makeURL, getRootPath

def buildSubfamilesFromGFF(job, gffID, halID, args):
    readGFFJob = Job.wrapJobFn(readGFF, halID=halID, gffID=gffID, args=args)

    buildSubfamiliesJob = Job.wrapJobFn(buildSubfamilies, elements=readGFFJob.rv(), args=args)
    writeGFFJob = Job.wrapJobFn(writeGFF, elements=buildSubfamiliesJob.rv())
    readGFFJob.addFollowOn(buildSubfamiliesJob)
    buildSubfamiliesJob.addFollowOn(writeGFFJob)
    job.addChild(readGFFJob)
    return writeGFFJob.rv()

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("hal", type=str)
    parser.add_argument("inGFF", type=str)
    parser.add_argument("outGFF", type=str)

    parser.add_argument("--graphsDir", type=str, default=None, help="Folder to store intermediate results")

    addRepeatAnnotatorOptions(parser)
    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()

    with Toil(args) as toil:
        halID = toil.importFile(makeURL(args.hal))
        gffID = toil.importFile(makeURL(args.inGFF))
        args.substMatrixID = toil.importFile(makeURL(os.path.join(getRootPath(), "blosum80.mat")))

        outGFFID = toil.start(Job.wrapJobFn(buildSubfamilesFromGFF, gffID=gffID, halID=halID, args=args))
        toil.exportFile(outGFFID, makeURL(args.outGFF))
