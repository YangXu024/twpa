import gdspy
import numpy as np
import math
import datetime
starttime = datetime.datetime.now()


def crosshair(x0,y0,w,l_x,l_y,layer_align):
    """
    Generates a crosshair marker at the specified position
    INPUT: center position, width and lengths of the arms
    OUTPUT: crosshair element
    """
    rv=gdspy.Rectangle((x0-w/2,y0-l_y/2),(x0+w/2,y0+l_y/2),layer=layer_align)
    rh=gdspy.Rectangle((x0-l_x/2,y0-w/2),(x0+l_x/2,y0+w/2),layer=layer_align)
    crosshair=gdspy.boolean(rv,rh,'or',layer=layer_align)
    return crosshair

def markers_ebeam_litho(center_coords,chip_side_x,chip_side_y,factor,layer_align):
    """
    Generates "crosshair" markers at the corners of each chip for alignment of the dicing saw and dots following that fo alignment of ebeam
    INPUT: array with center position of the chips, chip dimensions
    OUTPUT: dicing markers elements
    """
    elements=[]

    w=5*factor # width of the cross hair on the corner of a chip
    l_x= 100*factor # length of the cross hair
    l_y= 100*factor
    
    w_d = 20*factor # width of ebeam alignment dots
#    w_g = 10 # width of grid
    d_x = 300*factor # period of ebeam alignment dots in x and y direction
    d_y = 300*factor
    
    for c in center_coords:
#        print(c)
#        print(type(c))
        rv=gdspy.Rectangle((c[0]-chip_side_x/2-w/2,c[1]+chip_side_y/2-l_y/2),(c[0]-chip_side_x/2+w/2,c[1]+chip_side_y/2+l_y/2),layer=layer_align)
        rh=gdspy.Rectangle((c[0]-chip_side_x/2-l_x/2,c[1]+chip_side_y/2-w/2),(c[0]-chip_side_x/2+l_x/2,c[1]+chip_side_y/2+w/2),layer=layer_align)
        crosshair=gdspy.boolean(rv,rh,'or',layer=layer_align) #making a cross hair on the edge of a chip
        
        # rectangle for ebeam in vertical and horizontal direction
        reh=gdspy.Rectangle((c[0]-chip_side_x/2-w_d/2-(d_x+50),c[1]+chip_side_y/2-w_d/2),(c[0]-chip_side_x/2+w_d/2-(d_x+50),c[1]+chip_side_y/2+w_d/2),layer=layer_align)
        rev=gdspy.Rectangle((c[0]-chip_side_x/2-w_d/2,c[1]+chip_side_y/2+(d_y+50)-w_d/2),(c[0]-chip_side_x/2+w_d/2,c[1]+chip_side_y/2+(d_y+50)+w_d/2),layer=layer_align)
        
        # offsetting the horizontal dot to left and right
        offs_x=[0, -d_x,-2*d_x]
        for j in offs_x:
            elements.append(gdspy.copy(reh,j,0))
            
        offs_x=[2*d_x+100,3*d_x+100,4*d_x+100]
        for j in offs_x:
            elements.append(gdspy.copy(reh,j,0))
        # offsetting the vertical dot to top and bottom
        offs_y=[0, -2*d_y-100,-3*d_y-100,-4*d_y-100]
        for j in offs_y:
            elements.append(gdspy.copy(rev,0,j))
            
        offs_y=[d_y,2*d_y]
        for j in offs_y:
            elements.append(gdspy.copy(rev,0,j))
        
        offs_x=[0, chip_side_x]
        offs_y=[0, -chip_side_y]
        for j in offs_x:
            for k in offs_y:
                elements.append(gdspy.copy(crosshair,j,k))
    
    # array of 27 x 27 square dots of 10um side
    # square ebeam marker for chips number [0,6,last-6,last] last = len(center_coords)

#    positions = [3,9]
#    for j in positions:
#        se = gdspy.Rectangle((center_coords[j][0]+2500-w_g/2,center_coords[j][1]+6000-w_g/2),(center_coords[j][0]+2500+w_g/2,center_coords[j][1]+6000+w_g/2),layer=layer_align)
#        offs_x = [-1053,-966,-880,-795,-711,-628,-546,-465,-385,-306,-228,-151,-75,0,75,151,228,306,385,465,546,628,711,795,880,966,1053]
#        offs_y = [-1053,-966,-880,-795,-711,-628,-546,-465,-385,-306,-228,-151,-75,0,75,151,228,306,385,465,546,628,711,795,880,966,1053]
#        for i in offs_x:
#            for k in offs_y:
#                elements.append(gdspy.copy(se,i,k))
#    positions = [len(center_coords)-11,len(center_coords)-5]
#    for j in positions:
#        se = gdspy.Rectangle((center_coords[j][0]+2500-w_g/2,center_coords[j][1]-6000-w_g/2),(center_coords[j][0]+2500+w_g/2,center_coords[j][1]-6000+w_g/2),layer=layer_align)
#        offs_x = [-1053,-966,-880,-795,-711,-628,-546,-465,-385,-306,-228,-151,-75,0,75,151,228,306,385,465,546,628,711,795,880,966,1053]
#        offs_y = [-1053,-966,-880,-795,-711,-628,-546,-465,-385,-306,-228,-151,-75,0,75,151,228,306,385,465,546,628,711,795,880,966,1053]
#        for i in offs_x:
#            for k in offs_y:
#                elements.append(gdspy.copy(se,i,k))
    
    markers=gdspy.boolean(elements,None,'or',precision=0.001,max_points=199,layer=layer_align)
    return markers

def markers_ebeam_litho1(center_coords,chip_side_x,chip_side_y,layer_align):
    """
    Generates "crosshair" markers at the corners of each chip for alignment of the dicing saw and dots following that fo alignment of ebeam
    INPUT: array with center position of the chips, chip dimensions
    OUTPUT: dicing markers elements
    """
    elements=[]
    w=5 # width of the cross hair on the corner of a chip
    l_x= 100 # length of the cross hair
    l_y= 100
    
    for c in center_coords:
        print(c)
        print(type(c))
        rv=gdspy.Rectangle((c[0]-chip_side_x/2-w/2,c[1]+chip_side_y/2-l_y/2),(c[0]-chip_side_x/2+w/2,c[1]+chip_side_y/2+l_y/2),layer=layer_align)
        rh=gdspy.Rectangle((c[0]-chip_side_x/2-l_x/2,c[1]+chip_side_y/2-w/2),(c[0]-chip_side_x/2+l_x/2,c[1]+chip_side_y/2+w/2),layer=layer_align)
        crosshair=gdspy.boolean(rv,rh,'or',layer=layer_align) #making a cross hair on the edge of a chip
        elements.append(crosshair)
    return elements 
        
def markers_photo_litho(center_coords,chip_side_x,chip_side_y,w,l_x,l_y,layer_align):
    """
    Generates photo-litho markers that are cross hairs on the wafer edges with desired width and length
    INPUT: array with center position of the chips, chip dimensions
    OUTPUT: dicing markers elements
    """
    elements=[]
    # positions of photo-litho alignment markers
    positions=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
    for j in positions:
        elements.append(crosshair(center_coords[j][0],center_coords[j][1],w,l_x,l_y,layer_align))
        
    markers = gdspy.boolean(elements,None,'or',precision=0.001,max_points=199,layer=layer_align)
    return markers

##########################################################################################################################

def cap(inipoint, width, length, width_v, length_v, layer, theta, center=(0,0)):
    """
    Capacitor with two pins.

    Parameters
    ----------
    inipoint : initial point (x, y)
    width : width of capcitor
    length: length of capacitor
    layer: layer of capacitor

    Returns
    ----------
    polygon 

    Notes
    -----
    """
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length_v),
              (x0+width_v,y0+length_v),
              (x0+width_v,y0),
              (x0+width-width_v,y0),
              (x0+width-width_v,y0+length_v),
              (x0+width,y0+length_v),
              (x0+width,y0-length),
              (x0,y0-length),
              (x0,y0)]
    p = gdspy.Polygon(points,layer)
    return p.rotate(theta, center)
    
def cap_c(inipoint, width, length, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length),
              (x0+width,y0+length),
              (x0+width,y0),
              (x0,y0)]
    p = gdspy.Polygon(points,layer)
    return p.rotate(theta, center)

def cap_r(inipoint, width, length, width_v, length_v, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length+length_v),
              (x0+width,y0+length+length_v),
              (x0+width,y0+length_v),
              (x0+width_v,y0+length_v),
              (x0+width_v,y0),
              (x0,y0)]
    p = gdspy.Polygon(points,layer)
    
    return p.rotate(theta, center)

def cap_turn(inipoint, width, length, width_v, length_v, R, alpha, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points1 = [(x0,y0),
              (x0,y0+length_v),
              (x0-width_v,y0+length_v),
              (x0-width_v,y0-length),
              (x0,y0-length),
              (x0,y0)]
    p1 = gdspy.Polygon(points1,layer)
    p1.rotate(theta, center)
    
    #plot the cap rotating part
    r_inipoint = (x0-width_v,y0-R)#center of round
    p2 = gdspy.Round(r_inipoint,
                     R,
                     inner_radius=R-length,
                     initial_angle=0.5 * np.pi,
                     final_angle=0.5 * np.pi+alpha,
                     tolerance=0.00005)
    p2.rotate(theta, center)
    
    rot_inipoint = (x0-width_v/2,y0-R-width_v/2/math.tan(alpha/2))
    p3_1 = gdspy.Polygon(points1,layer)
    p3 = p3_1.rotate(alpha,rot_inipoint)
    p3.rotate(theta,center)
    
    p_b1 = gdspy.boolean(p1,p2,"or")
    p = gdspy.boolean(p_b1,p3,"or",layer = layer)
    return p

def cap_ri(inipoint, width, length, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length),
              (x0+width,y0+length),
              (x0+width,y0),
              (x0,y0)]
    p = gdspy.Polygon(points,layer)
    return p.rotate(theta, center)

def cap_sig(inipoint, width, length, width_v, length_v, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length_v),
              (x0+width_v,y0+length_v),
              (x0+width_v,y0),
              (x0+width,y0),
              (x0+width,y0-length),
              (x0,y0-length),
              (x0,y0)]
    p = gdspy.Polygon(points,layer)
    
    return p.rotate(theta, center)

def vias(inipoint, width, length, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length),
              (x0+width,y0+length),
              (x0+width,y0),
              (x0,y0)]
    p = gdspy.Polygon(points,layer)

    return p.rotate(theta, center)

def vias_turn(inipoint, width, length, width_v, ce_gap, R, alpha, layer, theta, center=(0,0)):
    x0 = inipoint[0]-width_v/2+width/2
    y0 = inipoint[1]+ce_gap-width/2
    points = [(x0,y0),
              (x0,y0+length),
              (x0-width,y0+length),
              (x0-width,y0),
              (x0,y0)]
    rot_inipoint = (inipoint[0]-width_v/2,inipoint[1]-R-width_v/2/math.tan(alpha/2))
    p_1 = gdspy.Polygon(points,layer)
    p = p_1.rotate(alpha,rot_inipoint)
    return p.rotate(theta, center)

def JJ(inipoint, radius, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    p = gdspy.Round((x0,y0),radius,tolerance=0.00005,layer=layer)
    return p.rotate(theta, center)

def etch(inipoint, width, length, width_e, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length),
              (x0+width,y0+length),
              (x0+width,y0),
              (x0,y0),
              (x0+width_e,y0+width_e),    
              (x0+width-width_e,y0+width_e),
              (x0+width-width_e,y0+length-width_e),
              (x0+width_e,y0+length-width_e),
              (x0+width_e,y0+width_e),
              (x0,y0)]
    p = gdspy.Polygon(points,layer)
    return p.rotate(theta, center)

def ind(inipoint, length_1st, length_spa, length_end, seg_1, seg_2, path_wid, layer, theta, center=(0,0)):    
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0-length_1st-path_wid,y0),
              (x0-length_1st-path_wid,y0+seg_1+path_wid),
              (x0-length_1st-path_wid-length_spa,y0+seg_1+path_wid),
              (x0-length_1st-path_wid-length_spa,y0),
              (x0-length_1st-path_wid*3-length_spa*2,y0),
              (x0-length_1st-path_wid*3-length_spa*2,y0+seg_2+2*path_wid),
              (x0-length_1st-path_wid*3-length_spa*2+path_wid+length_end,y0+seg_2+2*path_wid),
              (x0-length_1st-path_wid*3-length_spa*2+path_wid+length_end,y0+seg_2+path_wid),
              (x0-length_1st-path_wid*3-length_spa*2+path_wid,y0+seg_2+path_wid),
              (x0-length_1st-path_wid*3-length_spa*2+path_wid,y0+path_wid),
              (x0-length_1st-path_wid*2-length_spa,y0+path_wid),
              (x0-length_1st-path_wid*2-length_spa,y0+path_wid*2+seg_1),
              (x0-length_1st,y0+path_wid*2+seg_1),
              (x0-length_1st,y0+path_wid),
              (x0,y0+path_wid),
              (x0,y0)]
    p = gdspy.Polygon(points, layer)
    return p.rotate(theta, center)

def ind_end(inipoint, width, length, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0-length,y0),
              (x0-length,y0+width),
              (x0,y0+width),
              (x0,y0)]
    p = gdspy.Polygon(points, layer)
    return p.rotate(theta, center)
          
def gnd_ind(inipoint, width, length, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length),
              (x0+width,y0+length),
              (x0+width,y0)]
    p = gdspy.Polygon(points, layer)

    return p.rotate(theta, center)     

def gnd_bot(inipoint, width, length, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length),
              (x0+width,y0+length),
              (x0+width,y0)]
    p = gdspy.Polygon(points, layer)
    return p.rotate(theta, center)     

def rect(inipoint, width, length, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length),
              (x0+width,y0+length),
              (x0+width,y0)]
    p = gdspy.Polygon(points, layer)
    return p.rotate(theta, center) 

def gnd_m(inipoint, width, length, width_e, layer, theta, center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length),
              (x0+width,y0+length),
              (x0+width,y0),
              (x0,y0),
              (x0+width_e,y0+width_e),    
              (x0+width-width_e,y0+width_e),
              (x0+width-width_e,y0+length-width_e),
              (x0+width_e,y0+length-width_e),
              (x0+width_e,y0+width_e),
              (x0,y0)]
    p = gdspy.Polygon(points,layer)
    return p.rotate(theta, center) 

def pad_bool(inipoint,width1,width2,length1,length2,layer,theta,center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length1/2),
              (x0+width1,y0+length2/2),
              (x0+width1+width2,y0+length2/2),
              (x0+width1+width2,y0-length2/2),
              (x0+width1,y0-length2/2),
              (x0,y0-length1/2),
              (x0,y0)]
    p = gdspy.Polygon(points, layer)
    return p.rotate(theta, center)  

def pad_sig(inipoint,width1,width2,width3,width4,length1,length2,length3,layer,theta,center=(0,0)):
    x0 = inipoint[0]
    y0 = inipoint[1]
    points = [(x0,y0),
              (x0,y0+length1/2),
              (x0+width1,y0+length2/2),
              (x0+width1+width2,y0+length2/2),
              (x0+width1+width2+width3,y0+length3/2),
              (x0+width1+width2+width3+width4,y0+length3/2),
              (x0+width1+width2+width3+width4,y0-length3/2),
              (x0+width1+width2+width3,y0-length3/2),
              (x0+width1+width2,y0-length2/2),
              (x0+width1,y0-length2/2),
              (x0,y0-length1/2),
              (x0,y0)]
    p = gdspy.Polygon(points, layer)
    return p.rotate(theta, center)
   
if __name__ == "__main__":
    # Examples
    gdspy.current_library = gdspy.GdsLibrary()
    lib = gdspy.GdsLibrary()

###############################################################################
###############################################################################

#    gdspy.current_library = gdspy.GdsLibrary()
#    lib = gdspy.GdsLibrary()
    ##layer
    #0 Si
    #0 SiO2
    #0 JJ
    #0 Vias
    #0 metal Cap
    
    ###Parameter---------------------------------------------------------------
    factor = 4 
    c_gap = 3*factor                       # gap between 2 capacitor
    cc_gap = 3.5*factor                    # gap between cc and y center of e1 etching 
    cr_gap = 0.7*factor                    # boundary gap between cr and cr_inside
    
    # Capacitor 1
    c1_width = 13*factor
    c1_length = 14*factor
    c1_width_v = 3*factor
    c1_length_v = 6*factor
    c1_layer = 14
    
    # Capacitor 2
    c2_width = 13*factor
    c2_length = 14*factor
    c2_width_v = 3*factor
    c2_length_v = 6*factor
    c2_layer = 14
    
    # Capacitor 3
    c3_width = 13*factor
    c3_length = 9*factor
    c3_width_v = 3*factor
    c3_length_v = 6*factor
    c3_layer = 14
    
    #etching1: cap1
    e1_width = 13*factor
    e1_length = 7*factor  
    e1_gap = 1*factor
    e1_layer = 11
    
    #etching2: cap2
    e2_width = 13*factor
    e2_length = 7*factor
    e2_gap = 1*factor
    e2_layer = 11
    
    #etching_c: large ethching area for coupling capacitor
    ec_width = 13*factor
    ec_length = 16*factor
    ec_gap = 1*factor # etching width
    ec_layer = 11
    
    # Coupling Capacitor
    cc_width = 9*factor
    cc_length = 6*factor
    cc_layer = 14
    
    # Resonance Capacitor
    cr_width = 14*factor
    cr_length = 63*factor
    cr_width_v = 6*factor
    cr_length_v = 5*factor
    cr_layer = 14
    
    # Resonance Capacitor plane inside
    cri_width = (cr_width - 2*cr_gap)
    cri_length= (cr_length - 2*cr_gap)
    cri_layer = 14
    
    # Resonance Inductor
    lr_seg_1 = 57*factor # length of first seg
    lr_seg_2 = 60*factor # length of first seg
    lr_width = 0.8*factor # path width
    lr_length_1st = 5*factor
    lr_length_spa = 3.2*factor # spacing of the inductance length
    lr_length_end = 9*factor
    lr_path_wid = 0.8*factor
    lr_layer = 18
    
    #Resonance inductor end (connect to the bot of cap)
    lrend_length = 5*factor
    lrend_layer = 17 # bot of cap
    
    #Resonance Inductor ground plane
    lrg_width = 15*factor
    lrg_length = 64*factor
    lrg_layer = 13   # 4 for ground
    
    #vias
    v_width = 2.4*factor
    v_length = 2.4*factor
    v_layer = 12
    v_gap = 1*factor    # the gap between the bottom of antenna (rec reaches out for vias)
                  # and the bot of vias
    #vias resonator
    vr_width = 1.2*factor
    vr_layer = 12
    
    #JJ
    JJ_radius = 0.5*factor
    JJ_gap = v_gap
    JJ_layer = 10
    ce_gap = (e1_length/2+v_gap) # y gap between inipoint of c1 and center of e1
    
    ###Position (width for x, length for y)------------------------------------
    ###The origin (0,0) is set as the center of etch center of cap2
    x0_lr = (c2_width/2+c_gap/2-cr_width_v/2+cr_gap)
    y0_lr = (cc_gap+cc_length+cr_length_v+cr_gap)
    c1_inipoint = (-c1_width-c_gap-c2_width/2, -ce_gap)
    c2_inipoint = (-c2_width/2, -ce_gap)
    c3_inipoint = (c2_width/2+c_gap, -ce_gap)
    c_tn_inipoint = (-c_gap-c2_width/2, -ce_gap)
    e1_inipoint = (-e1_width/2+c_gap/2+c3_width+c_gap+c2_width/2, -e1_length/2)
    e2_inipoint = (-e2_width/2-c_gap/2-c2_width/2, -e2_length/2)
    ec_inipoint = (-ec_width/2+c_gap/2+c2_width/2, -e1_length/2)
    cc_inipoint = (c2_width/2+c_gap/2-cc_width/2, cc_gap)
    cr_inipoint = (c2_width/2+c_gap/2-cr_width_v/2, cc_gap+cc_length)
    cri_inipoint = (c2_width/2+c_gap/2-cr_width_v/2+cr_gap, cc_gap+cc_length+cr_length_v + cr_gap)
    lr_inipoint = (x0_lr, y0_lr)
    lrend_inipoint = (x0_lr, y0_lr+lr_width+lr_seg_2)
    lrg_inipoint = (c2_width/2+c_gap/2-cr_width_v/2-lrg_width, cc_gap+cc_length+cr_length_v+cr_length/2-lrg_length/2)
    v1_inipoint = (-c2_width/2-c_gap-c1_width+c1_width_v/2-v_width/2, -v_width/2)
    v2_inipoint = (-c2_width/2+c2_width_v/2-v_width/2, -v_width/2)
    v3_inipoint = (c2_width/2+c_gap+c3_width_v/2-v_width/2, -v_width/2)
    JJ1_inipoint = (-c2_width/2-c_gap-c1_width_v/2, 0)
    JJ2_inipoint = (c2_width/2-c2_width_v/2, 0)
    JJ3_inipoint = (c2_width/2+c_gap+c3_width-c2_width_v/2, 0)
    vr_inipoint = (x0_lr-lr_length_1st-2*lr_length_spa-2*lr_width+lr_length_end-vr_width,y0_lr+lr_width+lr_seg_2+lr_width/2-vr_width/2)
    
    ###paint-------------------------------------------------------------------
    #number of meander line
    num_mnl = math.floor(9/2)
    
    #One JTWPA unit: "+x" orientation of resonance cap
    theta = -0.5 * np.pi
    c1 = cap(c1_inipoint, c1_width, c1_length, c1_width_v, c1_length_v, c1_layer, theta)
    c2 = cap(c2_inipoint, c2_width, c2_length, c2_width_v, c2_length_v, c2_layer, theta)
    c3 = cap(c3_inipoint, c3_width, c3_length, c3_width_v, c3_length_v, c3_layer, theta)
    e1 = etch(e1_inipoint, e1_width, e1_length, e1_gap, e1_layer, theta)
    e2 = etch(e2_inipoint, e2_width, e2_length, e2_gap, e2_layer, theta)
    ec = etch(ec_inipoint, ec_width, ec_length, ec_gap, ec_layer, theta)
    v1=vias(v1_inipoint, v_width, v_length, v_layer, theta)
    v2=vias(v2_inipoint, v_width, v_length, v_layer, theta)
    v3=vias(v3_inipoint, v_width, v_length, v_layer, theta)
    vr = vias(vr_inipoint, vr_width, vr_width, vr_layer, theta)
    JJ1=JJ(JJ1_inipoint, JJ_radius, JJ_layer, theta)
    JJ2=JJ(JJ2_inipoint, JJ_radius, JJ_layer, theta)
    JJ3=JJ(JJ3_inipoint, JJ_radius, JJ_layer, theta)
    cc = cap_c(cc_inipoint, cc_width, cc_length, cc_layer, theta)
    cr = cap_r(cr_inipoint, cr_width, cr_length, cr_width_v, cr_length_v, cr_layer, theta)
    cri = cap_ri(cri_inipoint, cri_width, cri_length, cr_layer, theta)
    
    lr = ind(lr_inipoint, lr_length_1st, lr_length_spa, lr_length_end, lr_seg_1, lr_seg_2, lr_path_wid, lr_layer, theta)
    lrend = ind_end(lrend_inipoint, lr_width, lrend_length, lrend_layer, theta)
    lrg = gnd_ind(lrg_inipoint, lrg_width, lrg_length, lrg_layer, theta)
    lrg_bool = gdspy.boolean(lrg, vr, "not")
    lrg_bool = gdspy.boolean(lrg_bool, lrend, "not", layer = lrg_layer)
    
    ###Add to cell 
    cell = lib.new_cell('unit')
    cell.add([c1,c2,c3,e1,e2,ec,cc,cr,cri,lr,lrend,lrg_bool,v1,v2,v3,vr,JJ1,JJ2,JJ3])
    
    ### Draw-------------------------------------------------------------------
    #  Line
    p_1 = []
    Jline = lib.new_cell('Jline')
    unit_gap = 3*factor
    spacing = (c1_width+c2_width+c3_width+c_gap*2+unit_gap)
    spacing_len = c1_length+e1_length/2+cc_gap+cc_length+cr_length_v+cr_length
    # input turn left
    num_rr = 10
    R_turn1 = spacing*num_rr/np.pi-1
    num_turn = num_rr/2
    
    # rot new cell with cap_turn ##
    cell_turn = lib.new_cell('unit_turn')
    alpha = np.pi/num_rr
    R = (c1_width-2*c1_width_v)/alpha + c1_length/2
    c_tn = cap_turn(c_tn_inipoint, c1_width, c1_length, c1_width_v, c1_length_v, R, alpha, c1_layer, theta, center=(0,0))
    # v1r, the rot vias in the unit 
    v1r_inipoint = (-c2_width/2-c_gap-c1_width_v+c1_width_v/2-v_width/2, -v_width/2)
    v1r = vias_turn(c_tn_inipoint, v_width, v_length, c1_width_v, ce_gap, R, alpha, v_layer, theta)
    cell_turn.add([c_tn,c2,c3,e1,e2,ec,cc,cr,cri,lr,lrend,lrg_bool,v1r,v2,v3,vr,JJ1,JJ2,JJ3])
    
    vt = 1.95*factor #vertical comp
    tg = c_gap + (c1_width_v-v_width)/2 + 0.3*factor #tangent comp
    dx = 0 
    dy = 0
    dxf = 0 
    dyf = 0
    for i in range(math.ceil(num_turn)):
        alpha_cor = (i+1)/num_rr*np.pi
        dxf = dxf-vt*np.sin(alpha_cor) + tg*np.cos(alpha_cor)
        dyf = dyf+vt*np.cos(alpha_cor) + tg*np.sin(alpha_cor)
    
    
    for i in range(math.ceil(num_turn)):
        alpha = i/num_rr*180
        alpha_r = i/num_rr*np.pi
        alpha_cor = (i+1)/num_rr*np.pi
        p = gdspy.CellArray(cell_turn,1,1,(0,0), origin=(0,0),rotation=-90+alpha)
        p_1.append(p.translate(-R_turn1+R_turn1*np.sin(alpha_r)+dx-dxf, -R_turn1-spacing+R_turn1*(1-np.cos(alpha_r))+dy-dyf+spacing))
        dx = dx-vt*np.sin(alpha_cor) + tg*np.cos(alpha_cor)
        dy = dy+vt*np.cos(alpha_cor) + tg*np.sin(alpha_cor)
        
    # input

    num_in = 8
    x0_line = -spacing*num_in-R_turn1-dx # beginning of line (mid of first unit)
    y0_line = -R_turn1-dy
    inputline = gdspy.CellArray(cell,1,num_in,(0,spacing),origin=(x0_line, y0_line),rotation = -90)
    p_1.append(inputline)
    
    # line1
    num_line1 = 28
    num_line_bot = 40
    line1 = gdspy.CellArray(cell,1,num_line1,(0,spacing))
    p_1.append(line1)
    
    #### rotation1 on the top
    num_rr = 10
    vt1 = 1.95*factor #vertical comp
    tg1 =c_gap + (c1_width_v-v_width)/2 #tangent comp
    R_rot = spacing*num_rr/np.pi
    dx1f = 0
    dy1f = 0
    for i in range(num_rr):
        alpha_cor = (i+1)/num_rr*np.pi
        dx1f = dx1f+tg1*np.sin(alpha_cor) + vt1*np.cos(alpha_cor)
        dy1f = dy1f-vt1*np.sin(alpha_cor) + tg1*np.cos(alpha_cor)
    vt2 = 1.95*factor #vertical comp
    tg2 =c_gap + (c1_width_v-v_width)/2 #tangent comp
    dx2f = 0
    dy2f = 0
    for i in range(num_rr):
        alpha_cor = (i+1)/num_rr*np.pi   
        dx2f = dx2f+tg2*np.sin(alpha_cor) + vt2*np.cos(alpha_cor)
        dy2f = dy2f+vt2*np.sin(alpha_cor) - tg2*np.cos(alpha_cor)
 
    # line2 out
#    p_2 = []
    num_line2 = num_line1+num_line_bot
    for i in range(num_mnl):
        p_line = gdspy.CellArray(cell,1,num_line2,(0,spacing),origin=(dx1f*(i+1)+dx2f*i,spacing*num_line1+dy1f),x_reflection=True)
        
        p_1.append(p_line.translate(R_rot*(4*(i+1)-2),0))
    
    # line3 in
    num_line3 = num_line2
    for i in range(num_mnl):
        p_line = gdspy.CellArray(cell,1,num_line3,(0,spacing),origin=((dx1f+dx2f)*(i+1),spacing*(-num_line3+num_line1)+dy1f+dy2f))
        p_1.append(p_line.translate(R_rot*(4*i+4),0))
    
    #### rotation3 on the bot
    R_rot = spacing*num_rr/np.pi
    for j in range(num_mnl):
        dx2 = 0
        dy2 = 0
        for i in range(num_rr):
            alpha = i/num_rr*180
            alpha_r = i/num_rr*np.pi
            alpha_cor = (i+1)/num_rr*np.pi
            p = gdspy.CellArray(cell_turn,1,1,(0,0), origin=(0,0),rotation=180+alpha)
            p1 = p.translate(R_rot*(4*j+2)+R_rot*(1-np.cos(alpha_r))+dx2,(-num_line3+num_line1)*spacing-R_rot*np.sin(alpha_r)+dy2)
            p_1.append(p1.translate(dx1f*(j+1)+dx2f*j,dy1f))    
            dx2 = dx2+tg2*np.sin(alpha_cor) + vt2*np.cos(alpha_cor)
            dy2 = dy2+vt2*np.sin(alpha_cor) - tg2*np.cos(alpha_cor)
    
    #### rotation1 on the top continue
    for j in range(math.ceil(9/2)):
        dx1 = 0
        dy1 = 0
        for i in range(num_rr):
            alpha = i/num_rr*180
            alpha_r = i/num_rr*np.pi
            alpha_cor = (i+1)/num_rr*np.pi
            p = gdspy.CellArray(cell_turn,1,1,(0,0), origin=(0,0),rotation=180-alpha,x_reflection=True)
            p1 = p.translate(R_rot*4*j+R_rot*(1-np.cos(alpha_r))+dx1,num_line1*spacing+R_rot*np.sin(alpha_r)+dy1)
            p_1.append(p1.translate((dx1f+dx2f)*j,0))
            dx1 = dx1+tg1*np.sin(alpha_cor) + vt1*np.cos(alpha_cor)
            dy1 = dy1-vt1*np.sin(alpha_cor) + tg1*np.cos(alpha_cor)
    
    # line end 
    num_line_end = num_line1
    linelast = gdspy.CellArray(cell,1,num_line_end,(0,spacing),origin=(9*R_rot*2+dx1f*(num_mnl+1)+dx2f*num_mnl,num_line1*spacing+dy1f),x_reflection=True)
    p_1.append(linelast)
    
    # output turn 
    R_turn1 = spacing*num_rr/np.pi
    dx3f = 0
    dy3f = 0
    tg3 = c_gap + (c1_width_v-v_width)/2
    vt3 = 1.95*factor
    for i in range(math.ceil(num_turn)):
        alpha_cor = (i+1)/num_rr*np.pi
        dx3f = dx3f + tg3*np.sin(alpha_cor) + vt3*np.cos(alpha_cor)
        dy3f = dy3f + vt3*np.sin(alpha_cor) - tg3*np.cos(alpha_cor)
    dx3 = 0
    dy3 = 0
    for i in range(math.ceil(num_turn)):
        alpha = i/num_rr*180
        alpha_r = i/num_rr*np.pi
        alpha_cor = (i+1)/num_rr*np.pi
        p = gdspy.CellArray(cell_turn,1,1,(0,0), origin=(0,0),rotation=180+alpha)
        p_1.append(p.translate(9*R_rot*2+R_turn1*(1-np.cos(alpha_r))+dx3+dx1f*(num_mnl+1)+dx2f*num_mnl, -R_turn1*np.sin(alpha_r)+dy3+dy1f))
        dx3 = dx3 + tg3*np.sin(alpha_cor) + vt3*np.cos(alpha_cor)
        dy3 = dy3 + vt3*np.sin(alpha_cor) - tg3*np.cos(alpha_cor)
    
    # output
    num_out = 7
    xend_line = 9*R_rot*2+R_turn1+dx1f*(num_mnl+1)+dx2f*num_mnl+dx3f
    yend_line = -R_turn1+dy1f+dy3f
    lineout = gdspy.CellArray(cell,1,num_out,(0,spacing),origin=(xend_line, yend_line),rotation = -90)
    p_1.append(lineout)
    
    ### Ground-----------------------------------------------------------------
    #gnd bottom
    gb_inipoint = (0,0)
    gb_width = 8.5*factor
    gb_spac = 8.5*factor
    gb_layer = 15
    x0 = -spacing*num_in-R_turn1-spacing+gb_width*2
    y0 = -spacing*((num_line3-num_line1 + 3))-R_rot-5*gb_width*2
    gb = gnd_bot(gb_inipoint, gb_width, gb_width, gb_layer, theta)
    
    ### boolean for ground-----------------------------------------------------
    #start
    p1_b = []
    p2_b = []
    bool_width = spacing_len + 60*factor
    bool_metal_width = spacing_len +60*factor
    pcen_bias = spacing_len/2-c1_length-c1_length_v/2
    bool_inipoint = (x0,y0_line-pcen_bias)
    p_st = gdspy.Path(bool_width, bool_inipoint)
    p_st.segment(-R_turn1-x0,'+x')
    p_st.turn(R_turn1+pcen_bias,'l')
    p_st.segment((num_line1)*spacing,'+y')
    
    cell_gnd = lib.new_cell('cell_gnd')
    cell_gnd.add(gb)
    num_gb_wid = 235
    num_gb_len = 235
    gb_arr = gdspy.CellArray(cell_gnd, num_gb_wid, num_gb_len, (gb_width+gb_spac,gb_width+gb_spac), origin=(x0,y0))
    inv = gdspy.boolean(gb_arr, p_st, "not")
    Rbool_turn = R_turn1+pcen_bias+dx1f/2
    
    # middle
    for i in range(num_mnl):
        bool_inipoint1 = (-pcen_bias+R_rot*4*i+dx1f*i+dx2f*i, (num_line1-1)*spacing)
        p1 = gdspy.Path(bool_width, bool_inipoint1)
        p1.segment(spacing, '+y')
        p1.turn(Rbool_turn, 'rr')
        p1.segment((num_line3+0.5)*spacing, '-y')
        
        bool_inipoint2 = (-pcen_bias+(4*i+2)*R_rot+dx1f*(i+1)+dx2f*i, -(num_line3-num_line1-1)*spacing)
        p2 = gdspy.Path(bool_width, bool_inipoint2)
        p2.segment(spacing, '-y')
        p2.turn(Rbool_turn, 'll')
        p2.segment((num_line3)*spacing, '+y')
        
        inv = gdspy.boolean(inv, p1, "not")
        inv = gdspy.boolean(inv, p2, "not")
    
    # end
    bool_inipoint1 = (-pcen_bias+R_rot*4*(num_mnl)+dx1f*(num_mnl)+dx2f*num_mnl, (num_line1-1)*spacing)
    p1end = gdspy.Path(bool_width, bool_inipoint1)
    p1end.segment(spacing, '+y')
    p1end.turn(Rbool_turn, 'rr')
    p1end.segment((num_line1+0.5)*spacing, '-y')
    inv = gdspy.boolean(inv, p1end, "not")
    
    bool_inipoint2 = (-pcen_bias + 9*2*R_rot+dx1f*(num_mnl+1)+dx2f*num_mnl, spacing)
    p2end = gdspy.Path(bool_width, bool_inipoint2)
    p2end.segment(spacing, '-y')
    p2end.turn(Rbool_turn, 'l')
    p2end.segment((num_gb_wid*2-1)*gb_width+x0-R_rot*9*2-R_turn1, '+x')
    gndinv = gdspy.boolean(inv, p2end, "not",layer = gb_layer)
#    Jline.add(inv)
    
    ### Bool ground metal -----------------------------------------------------
    #gnd metal 
    gm_inipoint = (0,0)
    gm_spac = 3*factor
    gm_width = (gb_width+gb_spac)/3 - gm_spac
    gm_layer = 16
    gm_width_e = 0.3*factor
    gb_spacing_all = gb_width + gb_spac
    
    x0_gm1 = x0+gb_width+gb_spac/2-gm_width/2-gm_width-gm_spac 
    y0_gm1 = y0+gb_spac/2-gm_width/2 
    x0_gm2 = x0_gm1+gm_width+gm_spac 
    y0_gm2 = y0_gm1+gm_spac+gm_width
    
    x0_gm3 = x0_gm1+gm_width+gm_spac
    y0_gm3 = y0_gm2+gm_spac+gm_width
    num_gm_wid1 = num_gb_wid*3-3
    num_gm_len1 = num_gb_len-1
    num_gm_wid2 = num_gb_wid-1
    num_gm_len2 = num_gb_len-1
    num_gm_wid3 = num_gb_wid-1
    num_gm_len3 = num_gb_len-2
    inv = gdspy.boolean(gdspy.CellArray(cell_gnd, num_gb_wid, num_gb_len, (gb_width+gb_spac,gb_width+gb_spac), origin=(x0,y0)), p, "not")
    gm1 = gnd_m((x0_gm1,y0_gm1), gm_width, gm_width, gm_width_e, gm_layer, theta=0)
    gm2 = gnd_m((x0_gm2,y0_gm2), gm_width, gm_width, gm_width_e, gm_layer, theta=0)
    gm3 = gnd_m((x0_gm3,y0_gm3), gm_width, gm_width, gm_width_e, gm_layer, theta=0)
    
    stepx1 = gm_width+gm_spac
    stepy1 = gb_width+gb_spac
    stepx2 = gb_width+gb_spac
    stepy2 = gb_width+gb_spac
    stepx3 = gb_width+gb_spac
    stepy3 = gb_width+gb_spac
    
#    for i in range(num_mnl):
#        bool_inipoint1 = (-pcen_bias+R_rot*4*i+dx1f*i+dx2f*i, (num_line1-1)*spacing)
#        p1 = gdspy.Path(bool_metal_width, bool_inipoint1)
#        p1.segment(spacing, '+y')
#        p1.turn(Rbool_turn, 'rr')
#        p1.segment((num_line3+0.5)*spacing, '-y')
#        bool_inipoint2 = (-pcen_bias+(4*i+2)*R_rot+dx1f*(i+1)+dx2f*i, -(num_line3-num_line1-1)*spacing)
#        p2 = gdspy.Path(bool_metal_width, bool_inipoint2)
#        p2.segment(spacing, '-y')
#        p2.turn(Rbool_turn, 'll')
#        p2.segment((num_line3)*spacing, '+y')
#        p_st = gdspy.boolean(p_st, p1, "or")
#        p_st = gdspy.boolean(p_st, p2, "or")
#    p_st = gdspy.boolean(p_st, p1end, "or")
#    p_st = gdspy.boolean(p_st, p2end, "or")
#
#    for i in range(num_gm_wid1):
#        for j in range(num_gm_len1):
#            points = [(x0_gm1+i*stepx1,y0_gm1+j*stepy1)]       
#            if not gdspy.inside(points,p_st)[0]:
#                Jline.add(gdspy.copy(gm1,stepx1*i,stepy1*j))
#
#    for i in range(num_gm_wid2):
#        for j in range(num_gm_len2):
#            points1 = [(x0_gm2+i*stepx2,y0_gm2+j*stepy2)]
#            points2 = [(x0_gm3+i*stepx3,y0_gm3+j*stepy3)]
##            print(gdspy.inside(points1,p_st))
#            if not gdspy.inside(points1,p_st)[0]:
#                Jline.add(gdspy.copy(gm2,stepx2*i,stepy2*j))
#            if not gdspy.inside(points2,p_st)[0]:
#                Jline.add(gdspy.copy(gm3,stepx3*i,stepy3*j))

    ### Ground Padding ---------------------------------------------------------------
    # Padding horizontal
    ph_width = 1000*factor
    ph_gap = 1000*factor # distance from little squre to the whole bottom
    gnd_bound = 10*factor
    pad_layer = 14
    ph_length = (num_gb_len*2-1)*gb_width + ph_gap*2 + gnd_bound*2
    x0_p = x0-gnd_bound-ph_width # x0_p, y0_p left down corner point of the chip
    y0_p = y0-gnd_bound-ph_gap
    ph_inipoint = (x0_p,y0_p)
    ph = rect(ph_inipoint, ph_width, ph_length, layer=pad_layer, theta=0)
    phr_inipoint = (x0+gnd_bound+(num_gb_wid*2-1)*gb_width,y0_p)
    phr = rect(phr_inipoint, ph_width, ph_length, layer=pad_layer, theta=0)
    
    # Padding vertical
    pv_width = 2*gnd_bound+(num_gb_wid*2-1)*gb_width
    pv_length = ph_gap
    pv_inipoint = (x0-gnd_bound,y0_p)
    pv = rect(pv_inipoint, pv_width, pv_length-10*factor, layer=pad_layer, theta=0)
    pvu_inipoint = (x0-gnd_bound,y0_p+2*gnd_bound+(num_gb_len*2-1)*gb_width+pv_length)
    pvu = rect(pvu_inipoint, pv_width, pv_length, layer=pad_layer, theta=0)

    # bool pad input
    y0_sigIn = y0_line
    pin_wid = 20*factor
    pin_len = 20*factor
    pin_gap = 10*factor
    pin_con = pin_len + 2*pin_gap
    pin_len_s = 800*factor
    pin_gap_s = pin_len_s*pin_gap/pin_len
    pin_s = pin_len_s + pin_gap_s*2
    pbool_inipoint = (x0_p, y0_sigIn)
    p_bool = pad_bool(pbool_inipoint, ph_width-pin_wid,pin_wid,pin_s, pin_con, layer=pad_layer, theta=0)
    ph_inv = gdspy.boolean(ph, p_bool, "not", layer = pad_layer)
    
    #bool pad output
    y0_sigOut = y0_sigIn+spacing
    x_po_mir = x0_p+ph_width+gnd_bound+(num_gb_wid*2-1)*gb_width/2
    po_bool1 = p_bool.mirror((x_po_mir,1),(x_po_mir,2))
    po_bool = po_bool1.translate(0,yend_line-y0_line)
    phr_inv = gdspy.boolean(phr,po_bool,"not",layer = pad_layer)
    
    ### Signal Padding---------------------------------------------------------
    #add input signal pad
    cin_right = x0_line-c1_width/2-c_gap*2-c3_width # x position right hand side of cin 
    pad_inipoint = (x0_p,y0_sigIn)
    p_width3 = 1/3*pin_wid
    sig_wid = cin_right-x0_p-ph_width
    sig_len = (c1_length_v-ce_gap)*2
    pin_sig = pad_sig(pad_inipoint,ph_width-pin_wid,pin_wid*2/3,pin_wid/3,sig_wid,pin_len_s,pin_len,sig_len,layer=pad_layer,theta=0)
    
    cin_inipoint = (x0_line-c1_width/2-c_gap*2-c3_width,y0_line+ce_gap)
    cin_width = c1_width/2
    cin_length = c1_length
    cin_width_v = c1_width_v
    cin_length_v = c1_length_v
    cin = cap_sig(cin_inipoint, cin_width, cin_length, cin_width_v, cin_length_v, layer=c1_layer, theta=np.pi, center=cin_inipoint)
    pin_sig_cin = gdspy.boolean(pin_sig,cin,"or",layer = pad_layer)
    
    #input vias
    vin_inipoint = (cin_right-cin_width_v/2-v_width/2, y0_line-v_width/2)
    vin = vias(vin_inipoint, v_width, v_length, layer=v_layer, theta=0)

    ### add output signal pad
    cout_left = xend_line+c2_width/2+(num_out-1)*spacing+c_gap*2+c1_width # x position right hand side of cin 
    sig_len_out = (c1_length_v-ce_gap)*2
    sig_wid_out = cin_right-x0_p-ph_width+10*factor
    pout_sig1 = pad_sig(pad_inipoint,ph_width-pin_wid,pin_wid*2/3,pin_wid/3,sig_wid_out,pin_len_s,pin_len,sig_len_out,layer=pad_layer,theta=0)
    pout_sig2 = pout_sig1.mirror((x_po_mir,1),(x_po_mir,2))
    pout_sig = pout_sig2.translate(0,yend_line-y0_line)
  
    #add output coupling cap
    cout_inipoint = (cout_left,yend_line)
    cout_width = c1_width/2
    cout_length = c1_length
    cout_width_v = c1_width_v
    cout_length_v = c1_length_v
    cout1 = cap_sig(cout_inipoint, cout_width, cout_length, cout_width_v, cout_length_v, layer=c1_layer, theta=0)
    cout = cout1.mirror((x_po_mir,yend_line+ce_gap/2),(x_po_mir+1,yend_line+ce_gap/2))
    pout_sig_cout = gdspy.boolean(pout_sig,cout,"or",layer = pad_layer)
    p_2 = [pv,pvu,ph_inv,phr_inv,pin_sig_cin,pout_sig_cout]

    #vias out
    vout_inipoint = (cout_left+cin_width_v/2-v_width/2, yend_line-v_width/2)
    vout = vias(vout_inipoint, v_width, v_length, layer=v_layer, theta=0)
#    Jline.add(vout)
 
    #add out coupling cap etch
    eout_inipoint = (cout_left-c_gap/2+e1_width/2, yend_line+e1_length/2)
    eout_width = e1_width
    eout_length = e1_length
    eout_gap = e1_gap
    eout = etch(eout_inipoint, eout_width, eout_length, eout_gap, layer=e1_layer, theta=np.pi, center=eout_inipoint)
    p_3 = [vin,vout]
#    Jline.add(eout)

#    cellmarker.add(gdspy.CellArray(Jline, 8, 5, (9000,-spy), origin=(-32953,300+spy*2)))
    
#    asml = []
#    lib1 = gdspy.GdsLibrary()
#    lib1.read_gds('ASML_5500-300.gds')
#    for k,v in lib1.cells.items():
##        print(lib1.cells[k])
#        Jline.add(lib1.cells[k])
##        print(str(k))
#        print(k,v)
    
#    main = lib1.new_cell('main')   
#    main.add(gdspy.CellArray(Jline,1,1,(0,0),(x0_p,y0_p)))
#    
#    asml = lib1.cells['L_V_40_360_CHOP_5']
#    asml1 = lib1.extract(asml)
#    Jline.add(asml)
    lib1 = gdspy.GdsLibrary()
    lib.read_gds('ASML_5500-300.gds')
    asml = lib.cells['5500_300']
    # Load a GDSII file into a new library
#    gdsii = gdspy.GdsLibrary(infile='ASML_5500-300.gds')
    Jline.add(p_1)
    Jline.add(p_2)
    Jline.add(p_3)
    Jline.add(gndinv)
    Jline.add(eout)
    
# Use a previously-created library to load the file contents into
    chipx = 6000*factor
    chipy = 6000*factor
    chip_gap = 4000
    col = 3
    row = 3
    maskx = chipx*col+chip_gap*2
    masky = chipy*row+chip_gap*3
    mask_origin = (-x0_p-maskx/2,-y0_p-masky/2)
    maincell = lib.new_cell('main')
    maincell.add(gdspy.CellArray(Jline, col, row, (chipx+chip_gap,chipy+chip_gap), origin=mask_origin))
    maincell.add(asml)
#    lib.read_gds('ASML_5500-300.gds', units='convert')

### marker
    chip_side_x = 0
    chip_side_y = 0
    layer_align = 80
    cross_pty = 43000
    center_coords = [(-chipx-chip_gap,cross_pty)]
    cellmarker = lib.new_cell('marker')
    maincell.add(markers_ebeam_litho(center_coords,chip_side_x,chip_side_y,factor,layer_align))
 
#    gdspy.GdsLibrary.top_level()
#    lib.read_gds('filename.gds')
    lib.write_gds('Jline_withmarker_3.gds')