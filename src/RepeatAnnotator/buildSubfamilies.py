import argparse
import os
from toil.common import Toil
from toil.job import Job

from RepeatAnnotator.repeatAnnotator import addRepeatAnnotatorOptions, readGFF, writeGFF, buildSubfamiliesHeaviestBundling, buildSubfamiliesPartialOrderTreeBuilding, makeURL, getRootPath

def buildSubfamilesFromGFF(job, gffID, halID, args):
    readGFFJob = Job.wrapJobFn(readGFF, gffID = gffID)

    if args.heaviestBundling:
        buildSubfamiliesJob = Job.wrapJobFn(buildSubfamiliesHeaviestBundling, elements=readGFFJob.rv(), halID=halID, args=args)
    elif args.partialOrderTreeBuilding:
        buildSubfamiliesJob = Job.wrapJobFn(buildSubfamiliesPartialOrderTreeBuilding, elements=readGFFJob.rv(), halID=halID, args=args)
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
    addRepeatAnnotatorOptions(parser)
    Job.Runner.addToilOptions(parser)

    args = parser.parse_args()

    with Toil(args) as toil:
        halID = toil.importFile(makeURL(args.hal))
        gffID = toil.importFile(makeURL(args.inGFF))
        args.substMatrixID = toil.importFile(makeURL(os.path.join(getRootPath(), "blosum80.mat")))

        outGFFID = toil.start(Job.wrapJobFn(buildSubfamilesFromGFF, gffID=gffID, halID=halID, args=args))
        toil.exportFile(outGFFID, makeURL(args.outGFF))
