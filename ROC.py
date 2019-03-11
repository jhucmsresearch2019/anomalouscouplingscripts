import ROOT
import style
from array import array
import pprint

def D0Minus():

    print "Processing: D0Minus"
    
    #input files
    f1 = ROOT.TFile("~/work/JHUGen/JHUGenerator/STD_1.root")
    t1 = f1.Get("tree")
    
    f2 = ROOT.TFile("~/work/JHUGen/JHUGenerator/A3_1.root")
    t2 = f2.Get("tree")
    
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
    #Make ROC
    x = array('d', [0]*1001)
    y = array('d', [0]*1001)

    for i in range(1, 1000+1):
        x[i] = h1.Integral(1,i)
        y[i] = h2.Integral(1,i)
    
    g1 = ROOT.TGraph(1001, x, y)
    g1.SetTitle('ROC1')
    g1.SetLineColor(1);

    g1.GetXaxis().SetLimits(0, 1)
    g1.GetHistogram().SetMinimum(0)  
    g1.GetHistogram().SetMaximum(1)

    return g1

def ML(root_file_name, color):
    
    print "Processing: " + root_file_name 

    f = ROOT.TFile(root_file_name)
    t = f.Get("dataset/TestTree")
    
    #-------------------#
    #Fill histograms
    h3 = ROOT.TH1F("h3", "sig", 1000, -1, 1)
    h3.SetDirectory(0)
    h4 = ROOT.TH1F("h3", "background", 1000, -1, 1)
    h4.SetDirectory(0)

    for entry in t:
        if t.classID == 1:
            h3.Fill(t.BDTG)
        elif t.classID == 0:
            h4.Fill(t.BDTG)
        else: assert False

    h3.Scale(1/h3.Integral())
    h4.Scale(1/h4.Integral())
    
    #------------------#
    #Make ROC
    a = array('d', [0]*1001)
    b = array('d', [0]*1001)
    
    for i in range(1, 1000+1):
        a[i] = h3.Integral(1,i)
        b[i] = h4.Integral(1,i)
    
    g = ROOT.TGraph(1001, b, a)
    g.SetTitle(root_file_name)
    g.SetLineColor(color);
    
    return g


c = ROOT.TCanvas()



g1 = D0Minus();
g1.Draw("AL")

###ADD TMVA OUTPUT FILES HERE###

g2 = ML("./vbf/std_vs_a3_HJJpz.root", 2);
g3 = ML("./vbf/std_vs_a3_2.root", 3);
g4 = ML("./vbf/std_vs_a3_HJJpz.root", 4);


g2.Draw("L")
g3.Draw("L")
g4.Draw("L")

       
save_as_file_name = "D0Minus_D0Minues&Angles_Angles_Angles&JHHpzTest"
for ext in "png", "eps", "root", "pdf":
    c.SaveAs("~//www/ml/"+save_as_file_name+"." +ext)




