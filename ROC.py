import ROOT
import style
from array import array

def ROC():
    
#input files
    f1 = ROOT.TFile("~/work/JHUGen/JHUGenerator/STD_1.root")
    t1 = f1.Get("tree")
    
    f2 = ROOT.TFile("~/work/JHUGen/JHUGenerator/A3_1.root")
    t2 = f2.Get("tree")

    f3 = ROOT.TFile("./vbf/std_vs_a3_2.root")
    t3 = f3.Get("dataset/TestTree")
    
    f4 = ROOT.TFile("./vbf/std_vs_a3_3.root")
    t4 = f3.Get("dataset/TestTree")

#---------------------#
#Fill h1    
    h1 = ROOT.TH1F("h1", "STD", 1000, 0, 1)
    h1.SetDirectory(0)

    for entry in t1:
        h1.Fill(t1.D0minus)
        
    h1.Scale(1/h1.Integral())
#-------------------#
#Fill h2    
    
    h2 = ROOT.TH1F("h2", "STD", 1000, 0, 1)
    h2.SetDirectory(0)
    
    for entry in t2:
        h2.Fill(t2.D0minus)
    
    h2.Scale(1/h2.Integral())
#-------------------#
#Make ROC curve for h1 and h2
    
    x = array('d', [0]*1001)
    y = array('d', [0]*1001)

    for i in range(1000+1):
        x[i] = h1.Integral(1,i)
        y[i] = h2.Integral(1,i)

    c = ROOT.TCanvas()
    g1 = ROOT.TGraph(1001, x, y)
    g1.SetTitle('ROC1')
    g1.SetLineColor(2);
    g1.Draw("AL")
    
#-------------------#
#Fill h3 and h4        

    h3 = ROOT.TH1F("h3", "sig", 1000, -1, 1)
    h3.SetDirectory(0)
    h4 = ROOT.TH1F("h3", "background", 1000, -1, 1)
    h4.SetDirectory(0)

    for entry in t3:
        if t3.classID == 1:
            h3.Fill(t3.BDTG)
        elif t3.classID == 0:
            h4.Fill(t3.BDTG)
        else: assert False

    h3.Scale(1/h3.Integral())
    h4.Scale(1/h4.Integral())
#------------------#
#Make ROC for h3 and h4
    
    a = array('d', [0]*1001)
    b = array('d', [0]*1001)
    
    for i in range(1000+1):
        a[i] = h3.Integral(1,i)
        b[i] = h4.Integral(1,i)
    
    g2 = ROOT.TGraph(1001, b, a)
    g2.SetTitle('ROC2')
    g2.SetLineColor(4);
    g2.Draw("L")
#--------------------#
#Fill h5 and h6

    h5 = ROOT.TH1F("h5", "sig", 1000, -1, 1)
    h5.SetDirectory(0)
    h6 = ROOT.TH1F("h6", "background", 1000, -1, 1)
    h6.SetDirectory(0)
    
    for entry in t4:
        if t4.classID == 1:
            h5.Fill(t4.BDTG)
        elif t4.classID == 0:
            h6.Fill(t4.BDTG)
        else: assert False

    h5.Scale(1/h5.Integral())
    h6.Scale(1/h6.Integral())

#------------------#
#Make ROC for h5 and h6    

    e = array('d', [0]*1001)
    f = array('d', [0]*1001)

    for i in range(1000+1):
        e[i] = h5.Integral(1,i)
        f[i] = h6.Integral(1,i)
   
    g3 = ROOT.TGraph(1001, f, e)
    g3.SetTitle('ROC3')
    g3.SetLineColor(ROOT.kGreen+3);
    g3.Draw("L")
            
    save_as_file_name = "ROC3"

    for ext in "png", "eps", "root", "pdf":
        c.SaveAs("~//www/ml/"+save_as_file_name+"." +ext)


ROC()

