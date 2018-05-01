import unittest
import os
import shutil
from sonLib.bioio import getTempDirectory
from toil.job import Job
from toil.common import Toil

from RepeatAnnotator.repeatAnnotator import makeURL, runPoa2, parseHeaviestBundles, runTreeBuilding, getRootPath

fastaTestFile = """
>seq1
ACTCGACCACCAGTCA
>seq2
ACTCGACTACCAGTCA
>seq3
GTTCACTCACTCCTAACTAACTC"""

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
        self.assertEqual(results, frozenset([frozenset(['seq1', 'seq2']), frozenset(['seq3'])]))

