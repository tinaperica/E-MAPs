#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

from Tkinter import *
from tkFileDialog import askopenfilename
import tkMessageBox
import controller
import os.path

class DatasetLoader(Frame):
    def __init__(self, master,controller):
        Frame.__init__(self,master)
        self.controller = controller
        self.grid()
        self.label = Label(self,text="Dataset :")
        self.label.grid(column=0,row=0,sticky='EW')
        self.dataset = StringVar()
        self.dataset.set("None Loaded")
        self.entry = Entry(self,textvariable = self.dataset,state=DISABLED)
        self.entry.grid(column=1,row=0,sticky='EW')
        self.loadbutton = Button(self,text=u"Load",command=self.OnButtonClick)
        self.loadbutton.grid(column=2,row=0)
        self.label2 = Label(self,text="Details :")
        self.label2.grid(column=0,row=1,sticky='EW')
        self.labelVariable = StringVar()
        self.label = Label(self,
                              anchor="w",textvariable=self.labelVariable)
        self.label.grid(column=1,row=1,columnspan=2,sticky='EW')
        self.label1 = Label(self,text="Output :")
        self.label1.grid(column=0,row=2,sticky='EW')
        self.o = StringVar()
        self.o.set("Not set")
        self.output_entry = Entry(self,textvariable = self.o)
        self.output_entry.grid(column=1,row=2,columnspan=1,sticky='EW')
        self.grid_columnconfigure(1,weight=1)
    
    def OnButtonClick(self):
        self.loadbutton.configure(state=DISABLED)
        temp = askopenfilename(parent=self, title='Open file')
        if temp :
            try :
                self.controller.read_emap(temp)
                self.labelVariable.set("%s Alleles; %s Interactions" % (len(controller.genes),len(controller.ints)))
                self.dataset.set(temp)
                self.master.hasDataset()
                self.dirname,self.filename = os.path.split(temp)
                self.fileBaseName,extension = os.path.splitext(self.filename)
                self.o.set("%s_gia" % self.fileBaseName)
            except Exception , err:
                print err
                tkMessageBox.showerror(message='Error reading file: \n%s'%err)
        self.loadbutton.configure(state=ACTIVE)
    
    def get_output_details(self) :
        output_file = os.path.join(self.dirname,self.o.get())        
        return output_file

class MultipleTesting(Frame):
    def __init__(self, master,controller):
        Frame.__init__(self,master)
        self.controller = controller
        self.grid()
        self.bonf_val = StringVar()
        self.bonf_val.set("0.05")
        self.bonferonni = Entry(self,textvariable=self.bonf_val)
        self.bonferonni.grid(column=0,row=0)
        self.bonf_on = IntVar()
        self.bonf_on.set(1)
        self.bon_button = Checkbutton(self,text="Bonferroni",variable=self.bonf_on)
        self.bon_button.grid(column=0,row=1)
        self.fdr_val = StringVar()
        self.fdr_val.set("0.1")
        self.fdr = Entry(self,textvariable=self.fdr_val)
        self.fdr.grid(column=1,row=0)
        self.fdr_on = IntVar()
        self.fdr_on.set(0)
        self.fdr_button = Checkbutton(self,text="FDR",variable=self.fdr_on)
        self.fdr_button.grid(column=1,row=1)
        self.p_val = StringVar()
        self.p_val.set("0.001")
        self.p = Entry(self,textvariable=self.p_val)
        self.p.grid(column=2,row=0)
        self.p_on = IntVar()
        self.p_on.set(0)
        self.p_button = Checkbutton(self,text="Uncorrected",variable=self.p_on)
        self.p_button.grid(column=2,row=1)
    
    def get_details(self) :
        res = {}
        res['bonf'] = bool(self.bonf_on.get())     
        res['bonf_cutoff'] = float(self.bonf_val.get())     
        res['fdr'] = bool(self.fdr_on.get())     
        res['fdr_cutoff'] = float(self.fdr_val.get())
        res['rawp'] = bool(self.p_on.get())     
        res['rawp_cutoff'] = float(self.p_val.get())
        return res
        
class ComplexLoader(Frame):
    def __init__(self, master,controller):
        Frame.__init__(self,master)
        self.controller = controller
        self.grid()
        self.scrollbara = Scrollbar(self, orient=VERTICAL)
        self.lista = Listbox(self,selectmode=EXTENDED,yscrollcommand=self.scrollbara.set,exportselection=0)
        self.scrollbara.config(command=self.lista.yview)
        self.lista.grid(column=0,row=0,sticky="NSEW")
        self.scrollbara.grid(column=1,row=0,sticky="NS")
        self.scrollbarb = Scrollbar(self, orient=VERTICAL)
        self.listb = Listbox(self,selectmode=EXTENDED,yscrollcommand=self.scrollbarb.set,exportselection=0)
        self.scrollbarb.config(command=self.listb.yview)
        self.listb.grid(column=2,row=0,sticky="NSEW")
        self.scrollbarb.grid(column=3,row=0,sticky="NS")
        self.loada = Button(self,text=u"Load Annotations",command=self.LoadAnnotationsA)
        self.loada.grid(column=0,row=1,columnspan=1,sticky='EW')
        self.loadb = Button(self,text=u"Load Annotations",command=self.LoadAnnotationsB)
        self.loadb.grid(column=2,row=1,columnspan=1,sticky='EW')
        self.grid_columnconfigure(0,weight=1)
        self.grid_columnconfigure(2,weight=1)
        self.grid_rowconfigure(0,weight=1)
        self.loada.configure(state=DISABLED)    
        self.loadb.configure(state=DISABLED)
    
    def activate_buttons(self) :
        self.loada.configure(state=ACTIVE)    
        self.loadb.configure(state=ACTIVE)    
        
    def LoadAnnotationsA(self):
        self.loada.configure(state=DISABLED)
        temp = askopenfilename(parent=self, title='Open file')
        if temp :
            try :
                self.controller.read_complexes_a(temp)
                self.lista.delete(0,self.lista.size())
                for i in self.controller.get_complexes_a() :
                    self.lista.insert(END,i)
                self.lista.selection_set(0,self.lista.size()-1)
                if self.lista.size() > 0 and self.listb.size() > 0:
                    self.master.hasAnnotations()
            except Exception , err:
                print err
                tkMessageBox.showerror(message='Error reading file: \n%s'%err)
        self.loada.configure(state=ACTIVE)
    
    def LoadAnnotationsB(self):
        self.loadb.configure(state=DISABLED)
        temp = askopenfilename(parent=self, title='Open file')
        if temp :
            try :
                self.controller.read_complexes_b(temp)
                self.listb.delete(0,self.listb.size())
                for i in self.controller.get_complexes_b() :
                    self.listb.insert(END,i)
                self.listb.selection_set(0,self.listb.size()-1)
                if self.lista.size() > 0 and self.listb.size() > 0:
                    self.master.hasAnnotations()
            except Exception , err:
                print err
                tkMessageBox.showerror(message='Error reading file: \n%s'%err)
        self.loadb.configure(state=ACTIVE)  
    
    def get_selection_a(self):
        items = self.lista.curselection()
        return [self.lista.get(int(item)) for item in items]
    
    def get_selection_b(self):
        items = self.listb.curselection()
        return [self.listb.get(int(item)) for item in items]    

class GiaOptions(Frame):
    def __init__(self, master,controller):
        Frame.__init__(self,master)
        self.controller = controller
        self.grid()
        self.min_val = StringVar()
        self.min_val.set("3")
        self.min = Entry(self,textvariable=self.min_val)
        self.min.grid(column=1,row=0,sticky='W')
        self.min_on = IntVar()
        self.min_on.set(1)
        self.min_button = Label(self,text="Minimum measured interactions : ")
        self.min_button.grid(column=0,row=0,sticky='W')
        self.abs_on = IntVar()
        self.abs_on.set(0)
        Label(self,text="Use absolute interactions :").grid(column=0,row=1,sticky="W")
        self.abs_button = Checkbutton(self,text="",variable=self.abs_on)
        self.abs_button.grid(column=1,row=1,sticky='EW')        
        modeframe = Frame(self)
        modeframe.grid(column=0,row=2,columnspan=2,sticky='EW')
        Label(modeframe,text="Link type :").pack(side='left')    
        self.mode= StringVar()
        self.mode.set("both")
        self.both = Radiobutton(modeframe, text="Two Sided", variable=self.mode, value="both")
        self.both.pack(side='left')
        self.less = Radiobutton(modeframe, text="Less than", variable=self.mode, value="lt")
        self.less.pack(side='left')
        self.greater = Radiobutton(modeframe, text="Greater than", variable=self.mode, value="gt")
        self.greater.pack(side='left')
                        
class gia_app(Tk):
    def __init__(self,parent,controller):
        Tk.__init__(self,parent)
        self.controller = controller
        self.parent = parent
        self.initialize()
        self.minsize(self.winfo_reqwidth(),self.winfo_reqheight())

    def initialize(self):
        menubar = Menu(self)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Exit", command=self.quit)
        menubar.add_cascade(label="File", menu=filemenu)
        self.config(menu=menubar)
        
        self.loader = DatasetLoader(self,self.controller)
        self.loader.pack(expand=0,fill=X)
        Label(self,bg="black",fg="white",text="Annotations : ").pack(expand=0,fill=X,pady=10)
        self.cloader = ComplexLoader(self,self.controller)
        self.cloader.pack(expand=1,fill=BOTH)
        Label(self,bg="black",fg="white",text="Multiple testing correction : ").pack(expand=0,fill=X,pady=10)
        self.mtc = MultipleTesting(self,self.controller)
        self.mtc.pack()
        Label(self,bg="black",fg="white",text="Options: ").pack(expand=0,fill=X,pady=10)
        self.opts = GiaOptions(self,self.controller)
        self.opts.pack()
        self.runbutton = Button(self,text="Run",state=DISABLED,command=self.RunExperiment)
        self.runbutton.pack(fill=X,pady=10)
        self.statusLabelText= StringVar()
        self.statusLabelText.set("Waiting for genetic interaction dataset")
        self.statusLabel = Label(self,bg="lightgrey",fg="black",textvariable=self.statusLabelText)
        self.statusLabel.pack(expand=0,fill=X)
        
    def hasDataset(self):
       self.cloader.activate_buttons()
       self.statusLabelText.set("Waiting for annotations")
    
    def hasAnnotations(self):
      self.runbutton.configure(state=ACTIVE) 
      self.statusLabelText.set("Ready to run")
      
    def RunExperiment(self):
      self.runbutton.configure(state=DISABLED) 
      self.statusLabelText.set("Performing analysis")
      self.update()
      self.controller.run_experiment(self.mtc.get_details(),\
        self.loader.get_output_details(),self.cloader.get_selection_a(),self.cloader.get_selection_b(),bool(self.opts.abs_on.get()),\
        int(self.opts.min_val.get()),self.opts.mode.get())
      self.statusLabelText.set("Analysis complete")
      self.runbutton.configure(state=ACTIVE) 
           
if __name__ == "__main__":
    controller = controller.Controller()
    app = gia_app(None,controller)
    app.title('Genetic Interaction Analysis')
    app.mainloop()
