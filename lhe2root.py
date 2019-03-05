#!/usr/bin/env python

import argparse, os

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("outputfile")
  parser.add_argument("inputfile", nargs="+")
  g = parser.add_mutually_exclusive_group(required=True)
  g.add_argument("--vbf", action="store_true")
  g.add_argument("--zh", action="store_true")
  g.add_argument("--wh", action="store_true")
  parser.add_argument("--use-flavor", action="store_true")
  args = parser.parse_args()

  if os.path.exists(args.outputfile): raise IOError(args.outputfile+" already exists")
  for _ in args.inputfile:
    if not os.path.exists(_): raise IOError(_+" doesn't exist")

from array import array
import itertools

import ROOT

from lhefile import LHEFile_JHUGenVBFVH, LHEFile_Hwithdecay
from mela import TVar

bad = False
try:
  newf = ROOT.TFile(args.outputfile, "RECREATE")
  t = ROOT.TTree("tree", "tree")

  branchnames_float = "costheta1", "costheta2", "Phi1", "costhetastar", "Phi", "HJJpz"
  if args.zh or args.wh:
    branchnames_float += ("mV", "mVstar")
  if args.vbf:
    branchnames_float += ("q2V1", "q2V2")
  branchnames_float += ("pg1", "pg4", "pg1g4", "D0minus", "DCP")

  branchnames_int = ()

  branches = {name: array("f", [0]) for name in branchnames_float}
  branches_int = {name: array("i", [0]) for name in branchnames_int}

  for name in branchnames_float:
    t.Branch(name, branches[name], name+"/F")
  for name in branchnames_int:
    t.Branch(name, branches_int[names], name+"/I")
  branches.update(branches_int)

  if args.zh:
    g4 = 0.144057
  if args.wh:
    g4 = 0.1236136
  if args.vbf:
    g4 = 0.297979

  for inputfile in args.inputfile:
    print inputfile
    with LHEFile_JHUGenVBFVH(inputfile, isgen=args.use_flavor, reusemela=True) as f:
      for i, event in enumerate(f):
        if i != 0 and i % 1000 == 0:
          print "Processed", i, "events"

        if args.zh:
          process = TVar.Had_ZH
        elif args.wh:
          process = TVar.Had_WH
        elif args.vbf:
          process = TVar.JJVBF
        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = 1
        branches["pg1"][0] = event.computeProdP()

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz4 = g4
        branches["pg4"][0] = event.computeProdP()

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = 1
        event.ghz4 = g4
        branches["pg1g4"][0] = event.computeProdP() - branches["pg1"][0] - branches["pg4"][0]

        branches["D0minus"][0] = branches["pg1"][0] / (branches["pg1"][0] + branches["pg4"][0])
        branches["DCP"][0] = branches["pg1g4"][0] / (2 * (branches["pg1"][0] * branches["pg4"][0]) ** 0.5)

        if args.zh or args.wh:
          branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVHAngles(process)
          branches["mV"][0] = sum((particle.second for particle in event.associated), ROOT.TLorentzVector()).M()
          branches["mVstar"][0] = sum((particle.second for particle in itertools.chain(event.daughters, event.associated)), ROOT.TLorentzVector()).M()
        elif args.vbf:
          branches["q2V1"][0], branches["q2V2"][0], branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVBFAngles()
          branches["HJJpz"][0] = (event.mothers[0].second + event.mothers[1].second).Pz()

        t.Fill()

      print "Processed", i+1, "events"

  newf.Write()
except:
  bad = True
  raise
finally:
  if bad:
    try:
      os.remove(args.outputfile)
    except:
      pass
