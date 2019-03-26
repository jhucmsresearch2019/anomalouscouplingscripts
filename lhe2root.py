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
  parser.add_argument("--CJLST", action="store_true")
  parser.add_argument("--reweight-to", choices="fa3-0.5")
  args = parser.parse_args()

  if os.path.exists(args.outputfile): raise IOError(args.outputfile+" already exists")
  for _ in args.inputfile:
    if not os.path.exists(_) and not args.CJLST: raise IOError(_+" doesn't exist")

from array import array
import itertools

import ROOT

from lhefile import LHEFile_JHUGenVBFVH, LHEFile_Hwithdecay
from mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

def tlvfromptetaphim(pt, eta, phi, m):
  result = ROOT.TLorentzVector()
  result.SetPtEtaPhiM(pt, eta, phi, m)
  return result

class Event(object):
  doneinit = False
  def __init__(self, mela, daughters, associated, mothers, isgen, weight=1):
    self.daughters = daughters
    self.associated = associated
    self.mothers = mothers
    self.isgen = isgen
    self.mela = mela
    self.mela.setInputEvent(daughters, associated, mothers, isgen)
    self.doneinit = True
    self.weight = weight

  def __getattr__(self, attr):
    return getattr(self.mela, attr)

  def __setattr__(self, attr, value):
    if self.doneinit:
      return setattr(self.mela, attr, value)
    return super(Event, self).__setattr__(attr, value)

class CJLSTFile_VBFVH(object):
  __melas = {}
  def __init__(self, filename, *melaargs, **kwargs):
    self.isgen = kwargs.pop("isgen", True)
    reusemela = kwargs.pop("reusemela", False)
    if kwargs: raise ValueError("Unknown kwargs: " + ", ".join(kwargs))
    self.filename = filename
    if reusemela and melaargs in self.__melas:
      self.mela = self.__melas[melaargs]
    else:
      self.__melas[melaargs] = self.mela = Mela(*melaargs)

    self.f = None

  def __enter__(self):
    self.f = ROOT.TFile.Open(self.filename)
    self.t = self.f.Get("ZZTree/candTree")
    return self

  def __exit__(self, *args, **kwargs):
    return self.f.Close()

  def __iter__(self):
    for i, entry in enumerate(self.t):
      try:
        if self.isgen:
          raise NotImplementedError
        else:
          if len(entry.JetPt) < 2: continue
          if args.reweight_to == "fa3-0.5":
            if args.vbf:
              weight = entry.p_Gen_VBF_SIG_ghv1_1_JHUGen + g4**2 * entry.p_Gen_VBF_SIG_ghv4_1_JHUGen - g4 * entry.p_Gen_VBF_SIG_ghv1_1_ghv4_1_JHUGen
          elif args.reweightto is None:
            weight = 1
          yield Event(
            self.mela,
            SimpleParticleCollection_t(
              SimpleParticle_t(id, tlvfromptetaphim(pt, eta, phi, 0))
                for pt, eta, phi, id in itertools.izip(entry.LepPt, entry.LepEta, entry.LepPhi, entry.LepLepId)
            ),
            SimpleParticleCollection_t(
              SimpleParticle_t(id, tlvfromptetaphim(pt, eta, phi, m))
                for pt, eta, phi, m, id in itertools.chain(
                  itertools.izip(entry.JetPt, entry.JetEta, entry.JetPhi, entry.JetMass, (0, 0)),
#                  itertools.izip(entry.ExtraLepPt, entry.ExtraLepEta, entry.ExtraLepPhi, itertools.repeat(0), entry.ExtraLepLepId),
                )
            ),
            None,
            False,
            weight=weight,
          )
      except GeneratorExit:
        raise
      except:
        self.t.Show()
        raise
      finally:
        try:
          self.mela.resetInputEvent()
        except:
          pass

bad = False
try:
  newf = ROOT.TFile(args.outputfile, "RECREATE")
  t = ROOT.TTree("tree", "tree")

  branchnames_float = "costheta1", "costheta2", "Phi1", "costhetastar", "Phi", "HJJpz"
  if args.zh or args.wh:
    branchnames_float += ("mV", "mVstar")
  if args.vbf:
    branchnames_float += ("q2V1", "q2V2")
  branchnames_float += (
    "pg1", "pg4", "pg1g4", "D0minus", "DCP", "DCP_old",
    "pxH",  "pyH",  "pzH",  "EH",
    "pxj1", "pyj1", "pzj1", "Ej1",
    "pxj2", "pyj2", "pzj2", "Ej2",
    "weight",
  )

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

  fileclass = CJLSTFile_VBFVH if args.CJLST else LHEFile_JHUGenVBFVH

  for inputfile in args.inputfile:
    print inputfile
    with fileclass(inputfile, isgen=args.use_flavor, reusemela=True) as f:
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

        if branches["pg1"][0] == branches["pg4"][0] == 0:
          for _ in event.daughters: print _.first,; _.second.Print()
          for _ in event.associated: print _.first,; _.second.Print()
          continue

        branches["D0minus"][0] = branches["pg1"][0] / (branches["pg1"][0] + branches["pg4"][0])
        branches["DCP"][0] = branches["pg1g4"][0] / (2 * (branches["pg1"][0] * branches["pg4"][0]) ** 0.5)
        branches["DCP_old"][0] = branches["pg1g4"][0] / (branches["pg1"][0] + branches["pg4"][0])

        if args.zh or args.wh:
          branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVHAngles(process)
          branches["mV"][0] = sum((particle.second for particle in event.associated), ROOT.TLorentzVector()).M()
          branches["mVstar"][0] = sum((particle.second for particle in itertools.chain(event.daughters, event.associated)), ROOT.TLorentzVector()).M()
        elif args.vbf:
          branches["q2V1"][0], branches["q2V2"][0], branches["costheta1"][0], branches["costheta2"][0], branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]= event.computeVBFAngles()
          branches["HJJpz"][0] = sum((particle.second for particle in itertools.chain(event.daughters, event.associated)), ROOT.TLorentzVector()).Pz()

        pH = sum((particle.second for particle in event.daughters), ROOT.TLorentzVector())
        branches["pxH"][0] = pH.Px()
        branches["pyH"][0] = pH.Py()
        branches["pzH"][0] = pH.Pz()
        branches["EH"][0] = pH.E()

        pj1 = event.associated[0].second
        branches["pxj1"][0] = pj1.Px()
        branches["pyj1"][0] = pj1.Py()
        branches["pzj1"][0] = pj1.Pz()
        branches["Ej1"][0] = pj1.E()

        pj2 = event.associated[1].second
        branches["pxj2"][0] = pj2.Px()
        branches["pyj2"][0] = pj2.Py()
        branches["pzj2"][0] = pj2.Pz()
        branches["Ej2"][0] = pj2.E()

        branches["weight"][0] = event.weight

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
