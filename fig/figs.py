from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches

condext = loadtxt("../temp/ext_conds.txt")
exampleext = loadtxt("../temp/ext_example.txt")
funvalsext = loadtxt("../temp/ext_funs.txt")
potext = loadtxt("../temp/ext_pot.txt")
gradext = loadtxt("../temp/ext_grad.txt")
hessext = loadtxt("../temp/ext_hess.txt")
srcext = loadtxt("../temp/ext_src.txt")
targext = loadtxt("../temp/ext_targ.txt")

condint = loadtxt("../temp/int_conds.txt")
exampleint = loadtxt("../temp/int_example.txt")
funvalsint = loadtxt("../temp/int_funs.txt")
potint = loadtxt("../temp/int_pot.txt")
gradint = loadtxt("../temp/int_grad.txt")
hessint = loadtxt("../temp/int_hess.txt")
srcint = loadtxt("../temp/int_src.txt")
targint = loadtxt("../temp/int_targ.txt")

condextrtoz = loadtxt("../temp/extrtoz_conds.txt")
potextrtoz = loadtxt("../temp/extrtoz_pot.txt")
gradextrtoz = loadtxt("../temp/extrtoz_grad.txt")
hessextrtoz = loadtxt("../temp/extrtoz_hess.txt")

condintrtoz = loadtxt("../temp/intrtoz_conds.txt")
potintrtoz = loadtxt("../temp/intrtoz_pot.txt")
gradintrtoz = loadtxt("../temp/intrtoz_grad.txt")
hessintrtoz = loadtxt("../temp/intrtoz_hess.txt")

############################################
# COLORS, SIZES, STYLES
############################################

orange1 = (.667,.424,.224,1)
teal1 = (.113,.4,.4,1)
orange2 = (.667,.424,.224,.25)
teal2 = (.113,.4,.4,.25)
red1 = (.8,.1,.1,1)
clear = (0,0,0,0)
black1 = (0,0,0,1)
black2 = (0,0,0,.5)
grey1 = (0,0,0,.25)
figsize1 = (6,6)
figsize2 = (1200,400)

strn = '-'
stro = '-'
stre = ':'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# plot labels

newfunstrext = r'$\left (Q_n, K_n \right)$'
oldfunstrext = r'$\left (r^{-n}, K_n \right)$'
newfunstrint = r'$\left (r^n, P_n \right)$'
oldfunstrint = r'$\left (r^n, I_n \right)$'

newfunstrexta = r'$\left (Q_n, K_n \right), \lambda \to 0$'
oldfunstrexta = r'$\left (r^{-n}, K_n \right), \lambda \to 0$'
newfunstrinta = r'$\left (r^n, P_n \right), \lambda \to 0$'
oldfunstrinta = r'$\left (r^n, I_n \right), \lambda \to 0$'

newfunstrextr = r'$\left (Q_n, K_n \right), r \to 0$'
oldfunstrextr = r'$\left (r^{-n}, K_n \right), r \to 0$'
newfunstrintr = r'$\left (r^n, P_n \right), r \to 0$'
oldfunstrintr = r'$\left (r^n, I_n \right), r \to 0$'

# axis labels

lambdarstr = r'$\lambda R$'
kappastr = r'$\kappa $'
errpotstr = r'$E_u$'
errgradstr = r'$E_g$'
errhessstr = r'$E_h$'

# line widths for certain types of plots

condlw = 2
errlw = 3

# fontsizes for certain types of plots

condfs = 18
errfs = 24


##################################################
# PLOT GEOMETRY
##################################################

plt.rcParams.update({'font.size': condfs})

plt.close('all')

rad = 0.5
ts = linspace(0,2*pi,1000)
xcirc = rad*cos(ts)
ycirc = rad*sin(ts)

arrowscale=10

xt = targext[:,0]
yt = targext[:,1]
xs = srcext[:,0]
ys = srcext[:,1]
dx1 = srcext[:,4]/arrowscale
dy1 = srcext[:,5]/arrowscale
dx2 = srcext[:,6]/arrowscale
dy2 = srcext[:,7]/arrowscale

plt.figure(figsize=figsize1)
plt.plot(xt,yt,'X',markersize=10,color=orange2,markeredgecolor=orange1)
plt.plot(xs,ys,'o',markersize=8,color=teal2,markeredgecolor=teal1)
plt.plot(xcirc,ycirc,linewidth=2.0,color=black1)
plt.axis([-1,1,-1,1])
plt.axes().set_aspect('equal')
plt.tight_layout()
plt.savefig('geo_ext.pdf',format='pdf',bbox_inches='tight')

# newfig

plt.close('all')

xt = targint[:,0]
yt = targint[:,1]
xs = srcint[:,0]
ys = srcint[:,1]
dx1 = srcint[:,4]/arrowscale
dy1 = srcint[:,5]/arrowscale
dx2 = srcint[:,6]/arrowscale
dy2 = srcint[:,7]/arrowscale

plt.figure(figsize=figsize1)
plt.plot(xt,yt,'X',markersize=10,color=orange2,markeredgecolor=orange1)
plt.plot(xs,ys,'o',markersize=8,color=teal2,markeredgecolor=teal1)
plt.plot(xcirc,ycirc,linewidth=2.0,color=black1)
plt.axis([-1,1,-1,1])
plt.axes().set_aspect('equal')
plt.tight_layout()
plt.savefig('geo_int.pdf',format='pdf',bbox_inches='tight')

##################################################
# CONDITION NUMBER PLOTS
##################################################

plt.rcParams.update({'lines.linewidth': condlw})
plt.rcParams.update({'font.size': condfs})

plt.close('all')

condmbh = condext[:,2]
condnaive = condext[:,3]
condnaive[:] = [NaN if ele < 0 else ele for ele in condnaive]
condnaive[:] = [NaN if abs(ele) == 0.0 else ele for ele in condnaive] 
lam = condext[:,1]
scale = condext[:,4]
l = int(condext[0,0])
n = int(size(condext,0)/l)
condmbh = reshape(condmbh,(l,n),order='F')
condnaive = reshape(condnaive,(l,n),order='F')
lamsc = 0.5*reshape(multiply(lam,scale),(l,n),order='F')
lamsc = lamsc[0,:]
lamsc = lamsc[:]

ind = argsort(lamsc)

lamsc = lamsc[ind]

cmbh1 = condmbh[0,ind]
cmbh2 = condmbh[1,ind]
cmbh3 = condmbh[2,ind]
cmbh4 = condmbh[l-1,ind]

cnaive1 = condnaive[0,ind]
cnaive2 = condnaive[1,ind]
cnaive3 = condnaive[2,ind]
cnaive4 = condnaive[l-1,ind]

condmbhr = condextrtoz[:,2]
condnaiver = condextrtoz[:,3]
condnaiver[:] = [NaN if ele < 0 else ele for ele in condnaiver]
condnaiver[:] = [NaN if abs(ele) == 0.0 else ele for ele in condnaiver] 
lamr = condextrtoz[:,1]
scaler = condextrtoz[:,4]
lr = int(condextrtoz[0,0])
nr = int(size(condextrtoz,0)/lr)
condmbhr = reshape(condmbhr,(lr,nr),order='F')
condnaiver = reshape(condnaiver,(lr,nr),order='F')
lamscr = 0.5*reshape(multiply(lamr,scaler),(lr,nr),order='F')
lamscr = lamscr[0,:]
lamscr = lamscr[:]

indr = argsort(lamscr)

lamscr = lamscr[indr]

cmbh1r = condmbhr[0,indr]
cmbh2r = condmbhr[1,indr]
cmbh3r = condmbhr[2,indr]
cmbh4r = condmbhr[l-1,indr]

cnaive1r = condnaiver[0,indr]
cnaive2r = condnaiver[1,indr]
cnaive3r = condnaiver[2,indr]
cnaive4r = condnaiver[l-1,indr]

plt.figure(figsize=figsize1)
plt.plot(lamsc,cmbh1,color=teal1,label=newfunstrexta)
plt.plot(lamsc,cnaive1,color=orange1,label=oldfunstrexta)
plt.plot(lamsc,cmbh1r,'-x',color=teal1,label=newfunstrextr)
plt.plot(lamsc,cnaive1r,'-x',color=orange1,label=oldfunstrextr)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(kappastr)
plt.legend()
plt.savefig('cond_ext1.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,cmbh2,color=teal1,label=newfunstrexta)
plt.plot(lamsc,cnaive2,color=orange1,label=oldfunstrexta)
plt.plot(lamsc,cmbh2r,'-x',color=teal1,label=newfunstrextr)
plt.plot(lamsc,cnaive2r,'-x',color=orange1,label=oldfunstrextr)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(kappastr)

plt.legend()
plt.savefig('cond_ext2.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,cmbh3,color=teal1,label=newfunstrexta)
plt.plot(lamsc,cnaive3,color=orange1,label=oldfunstrexta)
plt.plot(lamsc,cmbh3r,'-x',color=teal1,label=newfunstrextr)
plt.plot(lamsc,cnaive3r,'-x',color=orange1,label=oldfunstrextr)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(kappastr)

plt.legend()
plt.savefig('cond_ext3.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,cmbh4,color=teal1,label=newfunstrexta)
plt.plot(lamsc,cnaive4,color=orange1,label=oldfunstrexta)
plt.plot(lamsc,cmbh4r,'-x',color=teal1,label=newfunstrextr)
plt.plot(lamsc,cnaive4r,'-x',color=orange1,label=oldfunstrextr)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(kappastr)

plt.legend()
plt.savefig('cond_ext4.pdf',format='pdf',bbox_inches='tight')

# interior

plt.close('all')

condmbh = condint[:,2]
condnaive = condint[:,3]
condnaive[:] = [NaN if ele < 0 else ele for ele in condnaive]
condnaive[:] = [NaN if abs(ele) == 0.0 else ele for ele in condnaive] 
lam = condint[:,1]
scale = condint[:,4]
l = int(condint[0,0])
n = int(size(condint,0)/l)
condmbh = reshape(condmbh,(l,n),order='F')
condnaive = reshape(condnaive,(l,n),order='F')
lamsc = 0.5*reshape(multiply(lam,scale),(l,n),order='F')
lamsc = lamsc[0,:]
lamsc = lamsc[:]

ind = argsort(lamsc)

lamsc = lamsc[ind]

cmbh1 = condmbh[0,ind]
cmbh2 = condmbh[1,ind]
cmbh3 = condmbh[2,ind]
cmbh4 = condmbh[l-1,ind]

cnaive1 = condnaive[0,ind]
cnaive2 = condnaive[1,ind]
cnaive3 = condnaive[2,ind]
cnaive4 = condnaive[l-1,ind]

condmbhr = condintrtoz[:,2]
condnaiver = condintrtoz[:,3]
condnaiver[:] = [NaN if ele < 0 else ele for ele in condnaiver]
condnaiver[:] = [NaN if abs(ele) == 0.0 else ele for ele in condnaiver] 
lamr = condintrtoz[:,1]
scaler = condintrtoz[:,4]
lr = int(condintrtoz[0,0])
nr = int(size(condintrtoz,0)/lr)
condmbhr = reshape(condmbhr,(lr,nr),order='F')
condnaiver = reshape(condnaiver,(lr,nr),order='F')
lamscr = 0.5*reshape(multiply(lamr,scaler),(lr,nr),order='F')
lamscr = lamscr[0,:]
lamscr = lamscr[:]

indr = argsort(lamscr)

lamscr = lamscr[indr]

cmbh1r = condmbhr[0,indr]
cmbh2r = condmbhr[1,indr]
cmbh3r = condmbhr[2,indr]
cmbh4r = condmbhr[l-1,indr]

cnaive1r = condnaiver[0,indr]
cnaive2r = condnaiver[1,indr]
cnaive3r = condnaiver[2,indr]
cnaive4r = condnaiver[l-1,indr]

plt.figure(figsize=figsize1)
plt.plot(lamsc,cmbh1,color=teal1,label=newfunstrinta)
plt.plot(lamsc,cnaive1,color=orange1,label=oldfunstrinta)
plt.plot(lamsc,cmbh1r,'-x',color=teal1,label=newfunstrintr)
plt.plot(lamsc,cnaive1r,'-x',color=orange1,label=oldfunstrintr)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(kappastr)

plt.legend()
plt.savefig('cond_int1.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,cmbh2,color=teal1,label=newfunstrinta)
plt.plot(lamsc,cnaive2,color=orange1,label=oldfunstrinta)
plt.plot(lamsc,cmbh2r,'-x',color=teal1,label=newfunstrintr)
plt.plot(lamsc,cnaive2r,'-x',color=orange1,label=oldfunstrintr)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(kappastr)

plt.legend()
plt.savefig('cond_int2.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,cmbh3,color=teal1,label=newfunstrinta)
plt.plot(lamsc,cnaive3,color=orange1,label=oldfunstrinta)
plt.plot(lamsc,cmbh3r,'-x',color=teal1,label=newfunstrintr)
plt.plot(lamsc,cnaive3r,'-x',color=orange1,label=oldfunstrintr)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(kappastr)

plt.legend()
plt.savefig('cond_int3.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,cmbh4,color=teal1,label=newfunstrinta)
plt.plot(lamsc,cnaive4,color=orange1,label=oldfunstrinta)
plt.plot(lamsc,cmbh4r,'-x',color=teal1,label=newfunstrintr)
plt.plot(lamsc,cnaive4r,'-x',color=orange1,label=oldfunstrintr)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(kappastr)

plt.legend()
plt.savefig('cond_int4.pdf',format='pdf',bbox_inches='tight')


##################################################
# ERROR PLOTS
##################################################

plt.rcParams.update({'lines.linewidth': errlw})
plt.rcParams.update({'font.size': errfs})

plt.close('all')

lam = potext[:,0]
scale = potext[:,5]
lamsc = 0.5*multiply(lam,scale)
ind = argsort(lamsc)
lamsc = lamsc[ind]

pot = potext[ind,1:5]
grad = gradext[ind,1:5]
hess = hessext[ind,1:5]

lamr = potextrtoz[:,0]
scaler = potextrtoz[:,5]
lamscr = 0.5*multiply(lamr,scaler)
indr = argsort(lamscr)
lamscr = lamscr[indr]

potr = potextrtoz[indr,1:5]
gradr = gradextrtoz[indr,1:5]
hessr = hessextrtoz[indr,1:5]

plt.figure(figsize=figsize1)
plt.plot(lamsc,pot[:,0],strn,color=teal1,label=newfunstrext)
plt.plot(lamsc,pot[:,1],stro,color=orange1,label=oldfunstrext)
plt.plot(lamsc,pot[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errpotstr)

plt.legend()
plt.savefig('errpot_ext.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)

plt.plot(lamscr,potr[:,0],strn,color=teal1,label=newfunstrext)
plt.plot(lamscr,potr[:,1],stro,color=orange1,label=oldfunstrext)
plt.plot(lamscr,potr[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errpotstr)

plt.legend()
plt.savefig('errpot_extr.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,grad[:,0],strn,color=teal1,label=newfunstrext)
plt.plot(lamsc,grad[:,1],stro,color=orange1,label=oldfunstrext)
plt.plot(lamsc,grad[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errgradstr)

plt.legend()
plt.savefig('errgrad_ext.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)

plt.plot(lamscr,gradr[:,0],strn,color=teal1,label=newfunstrext)
plt.plot(lamscr,gradr[:,1],stro,color=orange1,label=oldfunstrext)
plt.plot(lamscr,gradr[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errgradstr)

plt.legend()
plt.savefig('errgrad_extr.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,hess[:,0],strn,color=teal1,label=newfunstrext)
plt.plot(lamsc,hess[:,1],stro,color=orange1,label=oldfunstrext)
plt.plot(lamsc,hess[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errhessstr)

plt.legend()
plt.savefig('errhess_ext.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)

plt.plot(lamscr,hessr[:,0],strn,color=teal1,label=newfunstrext)
plt.plot(lamscr,hessr[:,1],stro,color=orange1,label=oldfunstrext)
plt.plot(lamscr,hessr[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errhessstr)

plt.legend()
plt.savefig('errhess_extr.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

# interior

plt.close('all')

lam = potint[:,0]
scale = potint[:,5]
lamsc = 0.5*multiply(lam,scale)
ind = argsort(lamsc)
lamsc = lamsc[ind]

pot = potint[ind,1:5]
grad = gradint[ind,1:5]
hess = hessint[ind,1:5]

lamr = potintrtoz[:,0]
scaler = potintrtoz[:,5]
lamscr = 0.5*multiply(lamr,scaler)
indr = argsort(lamscr)
lamscr = lamscr[indr]

potr = potintrtoz[indr,1:5]
gradr = gradintrtoz[indr,1:5]
hessr = hessintrtoz[indr,1:5]

plt.figure(figsize=figsize1)
plt.plot(lamsc,pot[:,0],strn,color=teal1,label=newfunstrint)
plt.plot(lamsc,pot[:,1],stro,color=orange1,label=oldfunstrint)
plt.plot(lamsc,pot[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errpotstr)

plt.legend()
plt.savefig('errpot_int.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)

plt.plot(lamscr,potr[:,0],strn,color=teal1,label=newfunstrint)
plt.plot(lamscr,potr[:,1],stro,color=orange1,label=oldfunstrint)
plt.plot(lamscr,potr[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errpotstr)

plt.legend()
plt.savefig('errpot_intr.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,grad[:,0],strn,color=teal1,label=newfunstrint)
plt.plot(lamsc,grad[:,1],stro,color=orange1,label=oldfunstrint)
plt.plot(lamsc,grad[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errgradstr)

plt.legend()
plt.savefig('errgrad_int.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)

plt.plot(lamscr,gradr[:,0],strn,color=teal1,label=newfunstrint)
plt.plot(lamscr,gradr[:,1],stro,color=orange1,label=oldfunstrint)
plt.plot(lamscr,gradr[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errgradstr)

plt.legend()
plt.savefig('errgrad_intr.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)
plt.plot(lamsc,hess[:,0],strn,color=teal1,label=newfunstrint)
plt.plot(lamsc,hess[:,1],stro,color=orange1,label=oldfunstrint)
plt.plot(lamsc,hess[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errhessstr)

plt.legend()
plt.savefig('errhess_int.pdf',format='pdf',bbox_inches='tight')

plt.close('all')

plt.figure(figsize=figsize1)

plt.plot(lamscr,hessr[:,0],strn,color=teal1,label=newfunstrint)
plt.plot(lamscr,hessr[:,1],stro,color=orange1,label=oldfunstrint)
plt.plot(lamscr,hessr[:,3],stre,color=orange1,label='Exact Difference')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlabel(lambdarstr)
plt.ylabel(errhessstr)

plt.legend()
plt.savefig('errhess_intr.pdf',format='pdf',bbox_inches='tight')

plt.close('all')


##################################################
# HEAT MAPS OF EXAMPLE FIELDS
##################################################

plt.rcParams.update({'font.size': errfs})

if (1 == 1):

    zheat1 = exampleext[:,0]
    zheat2 = exampleext[:,2]
    zheat3 = exampleext[:,4]
    xp = exampleext[:,8]
    yp = exampleext[:,9]
    mask = exampleext[:,6]
    mm = len(zheat1)
    m = int(sqrt(mm))
    zheat1 = reshape(zheat1,(m,m),order='F')
    zheat2 = reshape(zheat2,(m,m),order='F')
    zheat3 = reshape(zheat3,(m,m),order='F')
    xp = reshape(xp,(m,m),order='F')
    yp = reshape(yp,(m,m),order='F')
    yp = yp[:,1]
    xp = xp[1,:]
    mask = reshape(mask,(m,m),order='F')
    zheat1[mask > 100] = NaN
    zheat2[mask > 100] = NaN
    zheat3[mask > 100] = NaN
    
    plt.figure(figsize=figsize1)
    h = plt.pcolor(xp,yp,zheat1,vmin=nanmin(zheat1), vmax=nanmax(zheat1))
    h.cmap.set_under('white')
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.axes().set_aspect('equal')
    
    plt.savefig('examplepot_ext.png',format='png',bbox_inches='tight')
    
    plt.close('all')
    
    plt.figure(figsize=figsize1)
    h = plt.pcolor(xp,yp,zheat2,vmin=nanmin(zheat2), vmax=nanmax(zheat2))
    h.cmap.set_under('white')
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.axes().set_aspect('equal')
    plt.savefig('examplegx_ext.png',format='png',bbox_inches='tight')

    plt.close('all')

    plt.figure(figsize=figsize1)
    h = plt.pcolor(xp,yp,zheat3,vmin=nanmin(zheat3), vmax=nanmax(zheat3))
    h.cmap.set_under('white')
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.axes().set_aspect('equal')
    
    plt.savefig('examplehxy_ext.png',format='png',bbox_inches='tight')
    
    plt.close('all')
    
    # interior
    
    zheat1 = exampleint[:,0]
    zheat2 = exampleint[:,2]
    zheat3 = exampleint[:,4]
    xp = exampleint[:,8]
    yp = exampleint[:,9]
    mask = exampleint[:,6]
    mm = len(zheat1)
    m = int(sqrt(mm))
    zheat1 = reshape(zheat1,(m,m),order='F')
    zheat2 = reshape(zheat2,(m,m),order='F')
    zheat3 = reshape(zheat3,(m,m),order='F')
    xp = reshape(xp,(m,m),order='F')
    yp = reshape(yp,(m,m),order='F')
    yp = yp[:,1]
    xp = xp[1,:]
    mask = reshape(mask,(m,m),order='F')
    zheat1[mask > 100] = NaN
    zheat2[mask > 100] = NaN
    zheat3[mask > 100] = NaN
    
    plt.figure(figsize=figsize1)
    h = plt.pcolor(xp,yp,zheat1,vmin=nanmin(zheat1), vmax=nanmax(zheat1))
    h.cmap.set_under('white')
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.axes().set_aspect('equal')
    plt.savefig('examplepot_int.png',format='png',bbox_inches='tight')
    
    plt.close('all')
    
    plt.figure(figsize=figsize1)
    h = plt.pcolor(xp,yp,zheat2,vmin=nanmin(zheat2), vmax=nanmax(zheat2))
    h.cmap.set_under('white')
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.axes().set_aspect('equal')
    plt.savefig('examplegx_int.png',format='png',bbox_inches='tight')
    
    plt.close('all')
    
    plt.figure(figsize=figsize1)
    h = plt.pcolor(xp,yp,zheat3,vmin=nanmin(zheat3), vmax=nanmax(zheat3))
    h.cmap.set_under('white')
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.axes().set_aspect('equal')
    plt.savefig('examplehxy_int.png',format='png',bbox_inches='tight')

    plt.close('all')
