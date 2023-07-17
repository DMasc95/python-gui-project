import tkinter as tk
from tkinter.filedialog import askopenfilename
import os
import re
import numpy as np
from math import degrees

m = tk.Tk() #creates the main window
m.title('Brillouin Zone') #title of the window
m.geometry('600x490')
msg1 = tk.Label(m,text='Cell Paramaters')
msg1.place(x=60,y=10)
#msg2 = tk.Label(m,text='Atomic Species')
#msg2.place(x=60,y=40)
msg3 = tk.Label(m,text='Atomic Positions')
msg3.place(x=60,y=70)
#msg4 = tk.Label(m,text='K Points')
#msg4.place(x=60,y=100)


class Line:
    def __init__(self,slope,const):
        self.a = slope
        self.c = const

def findcross(L1,L2):
    global ponto
    x=(L1.c-L2.c)/(L2.a-L1.a)
    y=L1.a*x+L1.c
    ponto = (x,y)
        
valuek0x = tk.StringVar()
valuek0y = tk.StringVar()
valuenx = tk.StringVar()
valueny = tk.StringVar()
valuestep = tk.StringVar()
msgk0x = tk.Label(m,text="K0x")
msgk0x.place(x=60,y=360)
entry0x = tk.Entry(m,textvariable=valuek0x,width=8)
entry0x.place(x=90,y=360)
entry0x.bind('<Return>',(lambda event:Bz()))
msgk0y = tk.Label(m,text="K0y")
msgk0y.place(x=150,y=360)
entry0y = tk.Entry(m,textvariable=valuek0y,width=8)
entry0y.place(x=180,y=360)
entry0y.bind('<Return>',(lambda event: Bz()))
msgnx = tk.Label(m,text="Nx")
msgnx.place(x=60,y=400)
entry1 = tk.Entry(m,textvariable=valuenx,width=8)
entry1.place(x=90,y=400)
entry1.bind('<Return>',(lambda event:Bz()))
msgny = tk.Label(m,text="Ny")
msgny.place(x=150,y=400)
entry2 = tk.Entry(m,textvariable=valueny,width=8)
entry2.place(x=180,y=400)
entry2.bind('<Return>',(lambda event:Bz()))
msgstep = tk.Label(m,text="Step")
msgstep.place(x=240,y=400)
entry3 = tk.Entry(m,textvariable=valuestep,width=8)
entry3.place(x=280,y=400)
entry3.bind('<Return>',(lambda event:Bz()))




def parser(keyword, path):
    global posy,a1,a2,a3,b1,b2,b3,cell
    """
    This function receives as input a QE file and a keyword,
    and extracts the corresponding value
    """
    with open(path, "r") as f:
        content = f.readlines()

    strings = [
        "title",
        "calculation",
        "wf_collect",
        "outdir",
        "pseudo_dir",
        "prefix",
        "occupations",
        "smearing",
    ]
    numbers = ["ibrav", "nat", "ntyp", "nbnd"]
    fnumbers = ["ecutwfc", "degauss"]
    #IDEA: Make a parser with ply.
    #TODO: Future versions should raise appropriate warnings when the search fails.
    if keyword in strings:
        return re.search(keyword + r"\s*=\s*'(.+)'", content).group(1)
    if keyword in numbers:
        number = int(
            re.findall(r"\d+", os.popen("findstr " + keyword + " " + path).read())[0]
        )
        return number
    if keyword in fnumbers:
        number = float(
            re.findall(r"\d+\.\d+", os.popen("findstr " + keyword + " " + path).read())[
                0
            ]
        )
        return number
    if keyword == "CELL_PARAMETERS":
        for line in content:
            if line.find('alat') != -1:
                break
            if line.find('celldm(1)') != -1:
                linenr = content.index(line)
        celldm = int(re.findall(r"\d+",content[linenr])[-1])
        for line in content:
            if line.find(keyword) != -1:
                linenr = content.index(line)
        lines = (content[linenr+1].strip(),content[linenr+2].strip(),content[linenr+3].strip())
        cell = []
        for l in lines:
            vtrs = re.split(r"\s+",l)
            for i in vtrs:
                cell.append(float(i)*celldm)
        a1 = np.array(cell[0:3])
        a2 = np.array(cell[3:6])
        a3 = np.array(cell[6:9])
        b1 = 2*np.pi*((np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3))))
        b2 = 2*np.pi*((np.cross(a3,a1))/(np.dot(a1,np.cross(a2,a3))))
        b3 = 2*np.pi*((np.cross(a1,a2))/(np.dot(a1,np.cross(a2,a3))))
        tk.Label(m,text='b1 = ('+str(b1[0])+' , '+str(b1[1])+' , '+str(b1[2])+')').place(x=30,y=50)
        tk.Label(m,text='b2 = ('+str(b2[0])+' , '+str(b2[1])+' , '+str(b2[2])+')').place(x=30,y=90)
        tk.Label(m,text='b3 = ('+str(b3[0])+' , '+str(b3[1])+' , '+str(b3[2])+')').place(x=30,y=130)
        posy = 170
    elif keyword == "ATOMIC_SPECIES":
        #msg2.place(x=60,y=posy)
        for line in content:
            if line.find('ntyp') != -1:
                linenr = content.index(line)
        ntyp = int(re.findall(r"\d+",content[linenr])[0])
        for line in content:
            if line.find(keyword) != -1:    
                linenr = content.index(line)
        speci = []
        for i in range(1,ntyp+1):
            speci.append(content[linenr+i].strip())
        for i in range(0,ntyp):
            posy += 40
            tk.Label(m,text=str(speci[i])).place(x=30,y=posy)
    elif keyword == "ATOMIC_POSITIONS":
        posy += 40
        msg3.place(x=60,y=posy)
        for line in content:
            if line.find('nat') != -1:
                linenr = content.index(line)
        nat = int(re.findall(r"\d+",content[linenr])[0])
        for line in content:
            if line.find(keyword) != -1:
                linenr = content.index(line)
        atoms = []
        for i in range(1,nat+1):
            atoms.append(content[linenr+i].strip())
        for i in range(0,nat):
            posy += 40
            tk.Label(m,text=str(atoms[i])).place(x=30,y=posy)
    elif keyword == "K_POINTS":
        posy += 40
        #msg4.place(x=60,y=posy)
        for line in content:
            if line.find(keyword) != -1:
                linenr = content.index(line)
        kpoints = content[linenr+1]
        tk.Label(m,text=str(kpoints)).place(x=30,y=posy+40)       

def AboutWindow(): #definir a janela que fala sobre o programa
    a = tk.Toplevel(m)
    a.title('About')
    a.geometry('500x500')
    about_message = 'cenas' #arranjar maneira de meter aqui um texto de um ficheiro à parte, tipo pdf
    about = tk.Message(a,text=about_message)
    about.pack()

def HelpWindow():
    h = tk.Toplevel(m)
    h.title('Help')
    h.geometry('300x300')
    help_message = 'situações'
    help = tk.Message(h,text=help_message)
    help.pack()

def openfile():
    global path
    path = askopenfilename(filetypes=[('text files','*.txt'),('in files','*.in'),('all files','*.')])
    if os.path.exists(path):
        parser('CELL_PARAMETERS',path)
        #parser('ATOMIC_SPECIES',path)
        parser('ATOMIC_POSITIONS',path)
        #parser('K_POINTS',path)

def Bz(): #creates the Brilluoin Zone
    global x,b,y,slopes,blist
    bz = tk.Toplevel(m)
    bz.title('Brillouin Zone')
    bz.geometry('1000x500')
    bc = tk.Canvas(bz,height=2000,width=2000,bg='white')
    bc.pack()
    k0x = float(valuek0x.get())
    k0y = float(valuek0y.get())
    nx = int(valuenx.get())
    ny = int(valueny.get())
    step = float(valuestep.get())
    canvas_height = 500
    canvas_width = 1000
    cntrx = float(canvas_width/2)
    cntry = float(canvas_height/2)
    multiplier = 50
    b1[0] = b1[0]*multiplier
    b1[1] = b1[1]*multiplier
    b2[0] = b2[0]*multiplier
    b2[1] = b2[1]*multiplier
    for i in range(0,int(nx)):
        for l in range(0,int(ny)): 
            bc.create_oval(cntrx+k0x,cntry,cntrx+k0x,cntry,width=2,)
            bc.create_oval(cntrx+k0x,cntry-k0y,cntrx+k0x,cntry-k0y,width=2)
            k0y += step*multiplier
        k0x += step*multiplier
        k0y = float(valuek0y.get())
    bc.create_line(cntrx,cntry,cntrx+b1[0],cntry+b1[1],width=2,fill="green")
    bc.create_line(cntrx,cntry,cntrx+b2[0],cntry+b2[1],width=2,fill="green")
    theta = int(degrees(np.arccos(np.dot(b1,b2)/(np.linalg.norm(b1)*np.linalg.norm(b2)))))
    def draw(centerx,centery,extremos,width,color):
        bc.create_line(centerx+extremos[0][0],centery+extremos[0][1],centerx+extremos[1][0],centery+extremos[1][1],width=width,fill=color)
        bc.create_line(centerx+extremos[1][0],centery+extremos[1][1],centerx-extremos[0][0],centery+extremos[0][1],width=width,fill=color)
        bc.create_line(centerx-extremos[0][0],centery-extremos[0][1],centerx-extremos[1][0],centery-extremos[1][1],width=width,fill=color)
        bc.create_line(centerx-extremos[1][0],centery-extremos[1][1],centerx+extremos[0][0],centery-extremos[0][1],width=width,fill=color)
        bc.create_line(centerx+extremos[0][0],centery+extremos[0][1],centerx+extremos[0][0],centery-extremos[0][1],width=width,fill=color)
        bc.create_line(centerx-extremos[0][0],centery-extremos[0][1],centerx-extremos[0][0],centery+extremos[0][1],width=width,fill=color)
    def sqdraw(cntrx,cntry,M1,M2,width,color):
        bc.create_line(cntrx+M1[0],cntry+M2[1],cntrx+M1[0],cntry-M2[1],width=width,fill=color)
        bc.create_line(cntrx+M1[0],cntry+M2[1],cntrx-M1[0],cntry+M2[1],width=width,fill=color)
        bc.create_line(cntrx-M1[0],cntry-M2[1],cntrx-M1[0],cntry+M2[1],width=width,fill=color)
        bc.create_line(cntrx-M1[0],cntry-M2[1],cntrx+M1[0],cntry-M2[1],width=width,fill=color)
    if theta > 180:
        theta = 360-180
    if ((a1[0] != 0 and a1[1] == 0) or (a1[0] == 0 and a1[1] != 0)) and ((a2[0] != 0 and a2[1] == 0) or (a2[0] == 0 and a2[1] != 0)): 
        M1 = b1/2
        M2 = b2/2
        sqdraw(cntrx,cntry,M1,M2,2,"red")
        xtrvtrs = [b1,b2,-b1,-b2,b1+b2,b1-b2,-b1+b2,-b1-b2]
        for i in xtrvtrs:
            sqdraw(cntrx+i[0],cntry+i[1],M1,M2,2,"red")
    elif 0 < theta < 90:
        pontosm = [b1/2,b2/2,(-1/2)*(b1-b2),-b1/2,-b2/2,(1/2)*(b1-b2)]
        strpm = ['b1/2','b2/2','(-1/2)*(b1-b2)','-b1/2','-b2/2','(1/2)*(b1-b2)'] 
        xtrvtrs = [b1,b2,b2,-(b1-b2),-b1,-b2,(b1-b2)]
        slopes = []
        blist = []
        for i in strpm:
            x1 = cntrx
            y1 = cntry
            if ('-' and 'b1-b2') in i:
                x2 = cntrx-(b1[0]-b2[0])
                y2 = cntry-(b1[1]-b2[1])
            elif 'b1-b2' in i:
                x2 = cntrx+(b1[0]-b2[0])
                y2 = cntry+(b1[1]-b2[1])
            elif '-b1' in i:
                x2 = cntrx-b1[0]
                y2 = cntry-b1[1]
            elif 'b1' in i:
                x2 = cntrx+b1[0]
                y2 = cntry+b1[1]
            elif '-b2' in i:
                x2 = cntrx-b2[0]
                y2 = cntry-b2[1]
            elif 'b2' in i:
                x2 = cntrx+b2[0]
                y2 = cntry+b2[1]
            slope = (y2-y1)/(x2-x1) #declive dos vetores normais
            if slope == 0:
                slopes.append('inf')
                b = []
                blist.append(b)
            elif (x2-x1)==0:
                slopes.append(0.0)
                b=y2
                blist.append(b)
            else:
                invslope = -1/slope #declive dos vetores invertidos
                slopes.append(invslope) #acrescentar à lista dos declives das retas perpendiculares aos vetores
                nr = strpm.index(i)
                ind = pontosm[nr]
                b = ind[1] - ind[0]*invslope
                blist.append(b)
        print(slopes)
        print(blist)
        extremos = []
        for l in slopes:
            number = slopes.index(l)
            bnr = blist[number]
            if l == 0 or slopes[0] == 0:      #porque é que isto esta a dar 2 vezes o mesmo valor?????
                mpoint = pontosm[number]
                x = (mpoint[1]-blist[number+1])/slopes[number+1]
                extremos.append([x,mpoint[1]])
                continue
            elif l == 'inf' or slopes[0] == 'inf':
                mpoint = pontosm[number]
                y = (slopes[number-1])*(mpoint[0])+blist[number-1]
                extremos.append([mpoint[0],-y])
                continue
            elif number+1 == len(slopes):
                findcross(Line(l,bnr),Line(slopes[0],blist[0]))
                extremos.append([ponto[0],ponto[1]])
            elif slopes[number+1] == 0 or slopes[number+1] == 'inf':
                continue
            else:
                findcross(Line(l,bnr),Line(slopes[number+1],blist[number+1]))
                extremos.append([ponto[0],ponto[1]])
        draw(cntrx,cntry,extremos,2,"red")
        for i in xtrvtrs:
            draw(cntrx+i[0],cntry+i[1],extremos,2,"red")
    elif 90 < theta < 180:
        pontosm = [b1/2,(1/2)*(b1+b2),b2/2,-b1/2,(-1/2)*(b1+b2),-b2/2]
        strpm = ['b1/2','(1/2)*(b1+b2)','b2/2','-b1/2','(-1/2)*(b1+b2)','-b2/2'] 
        xtrvtrs = [b1,(b1+b2),b2,-b1,-(b1+b2),-b2]
        slopes = [] #lista dos declives das retas perpendiculares aos vetores
        blist = [] #lista dos bs das retas perpendiculares aos vetores
        for i in strpm:
            x1 = cntrx
            y1 = cntry
            if ('-' and 'b1+b2') in i:
                x2 = cntrx-(b1[0]+b2[0])
                y2 = cntry-(b1[1]+b2[1])
            elif 'b1+b2' in i:
                x2 = cntrx+(b1[0]+b2[0])
                y2 = cntry+(b1[1]+b2[1])
            elif '-b1' in i:
                x2 = cntrx-b1[0]
                y2 = cntry-b1[1]
            elif 'b1' in i:
                x2 = cntrx+b1[0]
                y2 = cntry+b1[1]
            elif '-b2' in i:
                x2 = cntrx-b2[0]
                y2 = cntry-b2[1]
            elif 'b2' in i:
                x2 = cntrx+b2[0]
                y2 = cntry+b2[1]
            slope = (y2-y1)/(x2-x1) #declive dos vetores normais
            if slope == 0:
                slopes.append('inf')
                b = []
                blist.append(b)
            elif (x2-x1)==0:
                slopes.append(0.0)
                b=y2
                blist.append(b)
            else:
                invslope = -1/slope #declive dos vetores invertidos
                slopes.append(invslope) #acrescentar à lista dos declives das retas perpendiculares aos vetores
                nr = strpm.index(i)
                ind = pontosm[nr]
                b = ind[1] - ind[0]*invslope
                blist.append(b)
        extremos = []
        for l in slopes:
            number = slopes.index(l)
            bnr = blist[number]
            if l == 0 or slopes[0] == 0:      #porque é que isto esta a dar 2 vezes o mesmo valor?????
                mpoint = pontosm[number]
                x = (mpoint[1]-blist[number+1])/slopes[number+1]
                extremos.append([x,mpoint[1]])
                continue
            elif l == 'inf' or slopes[0] == 'inf':
                mpoint = pontosm[number]
                y = (slopes[number-1])*(mpoint[0])+blist[number-1]
                extremos.append([mpoint[0],-y])
                continue
            elif number+1 == len(slopes):
                findcross(Line(l,bnr),Line(slopes[0],blist[0]))
                extremos.append([ponto[0],ponto[1]])
            elif slopes[number+1] == 0 or slopes[number+1] == 'inf':
                continue
            else:
                findcross(Line(l,bnr),Line(slopes[number+1],blist[number+1]))
                extremos.append([ponto[0],ponto[1]])
        draw(cntrx,cntry,extremos,2,"red")
        for i in xtrvtrs:
            draw(cntrx+i[0],cntry+i[1],extremos,2,"red")
    b1[0] = b1[0]/multiplier
    b1[1] = b1[1]/multiplier
    b2[0] = b2[0]/multiplier
    b2[1] = b2[1]/multiplier

btn1 = tk.Button(m,text='Run',width=15,command=Bz)
btn1.place(x=450,y=450)
btn1.bind('<Return>',(lambda event:Bz()))

menu = tk.Menu(m)
m.config(menu=menu)
menufile = tk.Menu(menu)
menu.add_cascade(label='File',menu=menufile)
menufile.add_command(label='Choose a file',command=openfile) 
menufile.add_separator()
menufile.add_command(label='Exit',command=m.quit)
help = tk.Menu(menu)
menu.add_cascade(label='Help',menu=help)
help.add_command(label='Help',command=HelpWindow)
help.add_command(label='About',command=AboutWindow)

m.mainloop()

#pwgui - ver software