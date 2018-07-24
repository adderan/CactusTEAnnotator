import unittest
import os
import shutil
from sonLib.bioio import getTempDirectory
from toil.job import Job
from toil.common import Toil

from RepeatAnnotator.repeatAnnotator import makeURL, runPoa, runTreeBuilding, getRootPath, Element, getDistances, joinByDistance

fastaTestFile = """
>seq1
ACTCGACCACCAGTCA
>seq2
ACTCGACTACCAGTCA
>seq3
GTTCACTCACTCCTAACTAACTC
>seq4
ACTCTCCCATTACACCATTCTACATAATCTTCAGGGCAGCAACATTCCAACACAATTACA"""

def treeBuildingTest(job, args):
    seqFile = job.fileStore.getLocalTempFile()
    with open(seqFile, 'w') as seqWrite:
        seqWrite.write(fastaTestFile)
    seqsID = job.fileStore.writeGlobalFile(seqFile)

    runPoaJob = Job.wrapJobFn(runPoa2, seqsID=seqsID, heaviestBundle=False, args=args)
    runTreeBuildingJob = Job.wrapJobFn(runTreeBuilding, graphID=runPoaJob.rv(), args=args)
    job.addChild(runPoaJob)
    runPoaJob.addFollowOn(runTreeBuildingJob)

    return runTreeBuildingJob.rv()

def heaviestBundlingTest(job, args):
    seqFile = job.fileStore.getLocalTempFile()
    with open(seqFile, 'w') as seqWrite:
        seqWrite.write(fastaTestFile)
    seqsID = job.fileStore.writeGlobalFile(seqFile)

    runPoaJob = Job.wrapJobFn(runPoa2, seqsID=seqsID, heaviestBundle=True, args=args)
    runParseHeaviestBundlesJob = Job.wrapJobFn(parseHeaviestBundles, runPoaJob.rv())
    
    job.addChild(runPoaJob)
    runPoaJob.addFollowOn(runParseHeaviestBundlesJob)
    return runParseHeaviestBundlesJob.rv()

def testPairwiseDistancesFn(job, args):
    elements = []

    seq1 = ">seq1\nACTACACATTACACACCACCAATAAATTTAAACACATTTACAC"
    seq1File = job.fileStore.getLocalTempFile()
    with open(seq1File, 'w') as seq1FileWrite:
        seq1FileWrite.write(seq1)
    elements.append(Element(chrom=None, start=None, end=None, group=0, name=1, seqID=job.fileStore.writeGlobalFile(seq1File)))

    seq2 = ">seq2\nACTACACATTACACACCACCAATAAATTTAAACACATTTACAC"
    seq2File = job.fileStore.getLocalTempFile()
    with open(seq2File, 'w') as seq2FileWrite:
        seq2FileWrite.write(seq2)
    elements.append(Element(chrom=None, start=None, end=None, group=0, name=2, seqID=job.fileStore.writeGlobalFile(seq2File)))
       
    seq3 = ">seq3\nACTACACATTACAAACGACCAATAAATTTCCACATATCTATCATTTACTACGTGTCGATCTAGAAACACATTTACAC"
    seq3File = job.fileStore.getLocalTempFile()
    with open(seq3File, 'w') as seq3FileWrite:
        seq3FileWrite.write(seq3)
    elements.append(Element(chrom=None, start=None, end=None, group=0, name=3, seqID=job.fileStore.writeGlobalFile(seq3File)))
    
    seq4 = ">seq4\nACTACACATTACAAATGACCAATAAATTTCCACATATCTATCATTTACTACGTGTCGATCTAGAAACACATTTACAC"
    seq4File = job.fileStore.getLocalTempFile()
    with open(seq4File, 'w') as seq4FileWrite:
        seq4FileWrite.write(seq3)
    elements.append(Element(chrom=None, start=None, end=None, group=0, name=4, seqID=job.fileStore.writeGlobalFile(seq4File)))

    args.distanceThreshold = 0.9

    pairwiseDistancesJob = Job.wrapJobFn(getDistances, elements=elements, args=args)
    joinByDistanceJob = Job.wrapJobFn(joinByDistance, elements=elements, distancesID=pairwiseDistancesJob.rv(), args=args)

    job.addChild(pairwiseDistancesJob)
    pairwiseDistancesJob.addFollowOn(joinByDistanceJob)

    return joinByDistanceJob.rv()

class RepeatAnnotatorTests(unittest.TestCase):
    def setUp(self):
        self.tempDir = getTempDirectory(os.getcwd())
        unittest.TestCase.setUp(self)
    def tearDown(self):
        shutil.rmtree(self.tempDir)
        unittest.TestCase.tearDown(self)

    def testHeaviestBundling(self):
        args = Job.Runner.getDefaultOptions(os.path.join(self.tempDir, "tmp_toil"))
        args.disableCaching = 'True'
        results = None
        with Toil(args) as toil:
            args.substMatrixID = toil.importFile(makeURL(os.path.join(getRootPath(), "blosum80.mat")))
            results = toil.start(Job.wrapJobFn(heaviestBundlingTest, args))
        print results
        self.assertEqual(results, frozenset([frozenset(['seq1', 'seq2']), frozenset(['seq3'])]))
    def testTreeBuilding(self):
        args = Job.Runner.getDefaultOptions(os.path.join(self.tempDir, "tmp_toil"))
        args.disableCaching = 'True'
        args.logLevel = "Debug"
        results = None
        with Toil(args) as toil:
            args.substMatrixID = toil.importFile(makeURL(os.path.join(getRootPath(), "blosum80.mat")))
            results = toil.start(Job.wrapJobFn(treeBuildingTest, args))
        self.assertEqual(results, frozenset([frozenset(['seq1', 'seq2']), frozenset(['seq3']), frozenset(['seq4'])]))
        print(results)

    def testPairwiseDistances(self):
        args = Job.Runner.getDefaultOptions(os.path.join(self.tempDir, "tmp_toil"))
        args.disableCaching = 'True'
        args.logLevel = "Debug"
        results = None

        with Toil(args) as toil:
            results = toil.start(Job.wrapJobFn(testPairwiseDistancesFn, args))
            self.assertEqual(len(results), 4)
            self.assertEqual(results[0].group, 0)
            self.assertEqual(results[1].group, 0)
            self.assertEqual(results[2].group, 1)
            self.assertEqual(results[3].group, 1)
