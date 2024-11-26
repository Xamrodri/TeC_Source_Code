################################################################################
# FUNCTIONS FILE
#
# 
# 2021/01/07
#
#
# Pascal Buri | High Mountain Glaciers and Hydrology | 
#  Swiss Federal Institute for Forest, Snow and Landscape Research, WSL |
#  Zürcherstrasse 111, 8903 Birmensdorf | pascal.buri@wsl.ch
#
#




################################################################################
# 'DISTRIBUTE_DAILYPRECIP': DISTRIBUTE DAILY PRECIP. SUMS TO SPECIFIC HOURS          
#
# DF_h:               data frame with hourly P values: 'Date' 'X' 'REF'
# DF_d:               data frame with daily P values: 'Date' 'X' 'REF'#   
# TARG_h:             vector with target hours (e.g. TARG_h<-18:23)
# X_ph:               hourly P for X (daily sums divided by 24)
# REF_ph:             hourly P for REF (daily sums divided by 24)

DISTRIBUTE_DAILYPRECIP<-function(DF_h,DF_d,TARG_h,X_ph,REF_ph){
  
  library(lubridate)
  
  idx<-which(format(DF_h$Date,'%Y-%m-%d') == as.character(DF_d[1,'Date']))
  idx2<-which(hour(DF_h[idx,'Date']) %in% TARG_h)
  
  DAT<-DF_h[idx,]
  DAT[-idx2,2]<-0
  DAT[-idx2,3]<-0
  DAT[idx2,2]<-X_ph
  DAT[idx2,3]<-REF_ph
  
  # DF_h[idx[-idx2],2]<-0
  # DF_h[idx[-idx2],3]<-0
  # DF_h[idx[idx2],2]<-X_ph
  # DF_h[idx[idx2],3]<-REF_ph
  # plot(DF_h[idx,2],type='l')
  # lines(DF_h[idx,3],col='darkgreen')
  rm(idx,idx2)
  
  return(DAT)
}





##########################################################################
# GENERATE SHAPEFILE FROM RASTER
#
#
# r:			Raster with e.g. 1 = area of interest & 0 = area outside
# projec:		Character string defining the projection, 
#				  e.g. for EPSG:3413 (WGS 84 / NSIDC Sea Ice Polar Stereographic North): 
#				  '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 
#				   +datum=WGS84 +units=m +no_defs'

RAST2POLY<-function(r,projec){
  
  xy<-coordinates(r)
  idx<-which(as.vector(r) > 0)
  xy<-xy[idx,]
  pts<-SpatialPoints(cbind(xy[,'x'],xy[,'y']))
  projection(pts)<-projec
  buf<-gBuffer(pts,byid=FALSE,width=res(r)[1])   
  outline<-gBuffer(buf,byid=FALSE,width=-res(r)[1]*0.5)
  outline<-unionSpatialPolygons(outline,rep(1,length(outline)))
  
  ##st_write(st_as_sf(outline),'outline.shp',driver="ESRI Shapefile",delete_dsn=T,quiet=T)
  
  return(outline)
}



##########################################################################
# DERIVE FLOW DIRECTION FROM SURFACE VELOCITY MAGNITUDE- 
#  & X-COMPONENT-RASTER
#
#  creates a raster with flow directions for each cell 
#  in geographical degree direction
#
# vv_r:            Raster of velocity magnitude
# vx_r:            Raster of velocity component in x-direction

FLOWDIR<-function(vv_r,vx_r){
  
  flowdir_r<-sin(vx_r/vv_r)*180/pi 
  #correct sheardirection raster values into geographic directions:
  #  currently: S=0?, E=90?, W=-90?, N=-180?/180?
  #  target:    S=180?, E=90?, W=270?, N=0?/360?)
  flowdir_v<-as.vector(flowdir_r)
  #W-sector
  idxW<-which(flowdir_v <= 0 & flowdir_v >= -180)
  flowdir_v[idxW]<-flowdir_v[idxW]*-1+180
  #E-sector
  idxE<-which(flowdir_v > 0 & flowdir_v <= 180)
  flowdir_v[idxE]<-180-flowdir_v[idxE]
  raster::values(flowdir_r)<-flowdir_v   #if values() not working: restart R
  return(flowdir_r)
}



##########################################################################
# CALCULATE DU & DV FROM RASTER (VELOCITY DIFFERENCE IN X & Y DIRECTION)
#
#  creates a list with four elements: dx, dy & two matrices du, dv 
#
# DX:              No. of neighbouring pixels to consider (1 or 2)
# vx_r:            Raster of velocity component in x-direction
# vy_r:            Raster of velocity component in y-direction

DUDV<-function(DX,vx_r,vy_r){
  
  #distance in x direction 
  dx<-res(vx_r)[1] * DX
  
  #distance in y direction
  dy<-dx
  
  ### Case 1: strain derived from 1 neighbouring cell ###
  if(DX == 1){
    #du  (velocity difference in x-direction)
    vx_m<-as.matrix(vx_r)
    dim_vx_m<-dim(vx_m)
    du_m<-matrix(nrow=dim_vx_m[1],ncol=dim_vx_m[2])
    for(ROW in 1:dim_vx_m[1]){
      du_v<-vector()
      for(COL in 1:(dim_vx_m[2]-1)){
        du_v[COL]<-vx_m[ROW,COL+1] - vx_m[ROW,COL]
      }
      du_m[ROW,]<-c(du_v,NA)   #add NA for raster margin
    }
    
    #dv  (velocity difference in y-direction)
    vy_m<-as.matrix(vy_r)
    dim_vy_m<-dim(vy_m)
    dv_m<-matrix(nrow=dim_vy_m[1],ncol=dim_vy_m[2])
    for(COL in 1:dim_vy_m[2]){
      dv_v<-vector()
      for(ROW in 2:dim_vy_m[1]){
        dv_v[ROW]<-vy_m[ROW-1,COL] - vy_m[ROW,COL]
      }
      dv_m[,COL]<-dv_v
    }
  }
  
  ### Case 2: strain derived from 2 neighbouring cells ###
  if(DX == 2){
    #du  (velocity difference in x-direction)
    vx_m<-as.matrix(vx_r)
    dim_vx_m<-dim(vx_m)
    du_m<-matrix(nrow=dim_vx_m[1],ncol=dim_vx_m[2])
    for(ROW in 1:dim_vx_m[1]){
      du_v<-vector()
      for(COL in 2:(dim_vx_m[2]-1)){
        dul<-(vx_m[ROW,COL] - vx_m[ROW,COL-1])  #left to center
        dur<-(vx_m[ROW,COL+1] - vx_m[ROW,COL])  #center to right
        du_v[COL]<-(dul+dur)/2                  #average
      }
      du_m[ROW,]<-c(du_v,NA)   #add NA for raster margin
    }
    
    
    #dv  (velocity difference in y-direction)
    vy_m<-as.matrix(vy_r)
    dim_vy_m<-dim(vy_m)
    dv_m<-matrix(nrow=dim_vy_m[1],ncol=dim_vy_m[2])
    for(COL in 1:dim_vy_m[2]){
      dv_v<-vector()
      for(ROW in 2:(dim_vy_m[1]-1)){
        dvb<-(vy_m[ROW,COL] - vy_m[ROW+1,COL])  #below to center
        dva<-(vy_m[ROW-1,COL] - vy_m[ROW,COL])  #center to above
        dv_v[ROW]<-(dvb+dva)/2                  #average
      }
      dv_m[,COL]<-c(dv_v,NA)   
    }
  }
  
  ls<-list(dx,dy,du_m,dv_m)
  names(ls)<-c('dx','dy','du','dv')
  return(ls)
}



##########################################################################
# CALCULATE PRINCIPAL STRAIN COMPONENTS
#
#
# du_m:           Velocity difference in x direction
# dv_m:           Velocity difference in y direction
# dx:             Distance in x direction 
# dy:             Distance in y direction
# r:              Raster to fill with values 

PRSTRAIN<-function(du_m,dv_m,dx,dy,r){
  
  #Exx (strain rate xx)
  Exx<-du_m/dx
  
  #Eyy (strain rate yy)
  Eyy<-dv_m/dy
  
  #Exy (strain rate xy)
  Exy<-0.5*(du_m/dy + dv_m/dx)
  
  ##derive min. and max. principal strain (Nye 1959, JG; Eqs. 4 & 5)
  lambda_min<-0.5*(Exx+Eyy) - sqrt(0.25*(Exx-Eyy)^2 + Exy^2)
  lambda_max<-0.5*(Exx+Eyy) + sqrt(0.25*(Exx-Eyy)^2 + Exy^2)
  
  lambda_min_r<-r
  raster::values(lambda_min_r)<-lambda_min
  names(lambda_min_r)<-'lambda_min'
  lambda_max_r<-r
  raster::values(lambda_max_r)<-lambda_max
  names(lambda_max_r)<-'lambda_max'
  
  #Principal strain direction (Nye 1959, JG; Eq. 6)
  # tan(2*Dir)= (2*Exy)/(Exx-Eyy)
  Dir<-atan((2*Exy)/(Exx-Eyy))/2
  Dir<-Dir*180/pi
  
  dir_r<-r
  raster::values(dir_r)<-Dir
  names(dir_r)<-'PrincipalDir'
  
  ls<-list(Exx,Eyy,Exy,lambda_min_r,lambda_max_r,dir_r)
  names(ls)<-c('Exx','Eyy','Exy','lambda_min','lambda_max','dir')
  return(ls)
}



##########################################################################
# EXTENSION/COMPRESSION VECTORS FROM PRINCIPAL STRAIN
#
#  creates dataframe with start/end-coordinates of 
#   extension and compression vectors
#   (needed for visualization)
#
# xy:              Matrix of xy-coordinates
# dir:             Vector of principal strain directions ([deg], rel. to x-axis)
# Lmx:             Vector of max. principal strain (lambda max.)
# Lmn:             Vector of min. principal strain (lambda min.)
# f:               Vector-enlargement factor (to plot vectors bigger)

ExtCom_vectors<-function(xy,dir,Lmx,Lmn,f){
  
  DF<-cbind.data.frame(xy,dir,Lmx,Lmn)
  DF$dir<-DF$dir*pi/180
  
  ### EXTENSION VECTOR COORDINATES ###
  #define prefactors for extension
  f_ext_sin<-sin(DF$dir)*DF$Lmx*f*0.5
  f_ext_cos<-cos(DF$dir)*DF$Lmx*f*0.5
  #Extension top point Y
  DF$Etop_y<-DF$y + f_ext_cos
  #Extension bottom point Y
  DF$Ebot_y<-DF$y - f_ext_cos
  idxpos<-which(DF$dir >= 0)
  idxneg<-which(DF$dir < 0)
  #Extension top point X
  DF[idxpos,'Etop_x']<-DF[idxpos,'x'] - f_ext_sin[idxpos]
  DF[idxneg,'Etop_x']<-DF[idxneg,'x'] + f_ext_sin[idxneg]
  #Extension bottom point X
  DF[idxpos,'Ebot_x']<-DF[idxpos,'x'] + f_ext_sin[idxpos]
  DF[idxneg,'Ebot_x']<-DF[idxneg,'x'] - f_ext_sin[idxneg]
  
  
  ### COMPRESSION VECTOR COORDINATES ###
  #define factor for compression
  f_com_sin<-sin(pi/2-DF$dir)*DF$Lmn*f*0.5
  f_com_cos<-cos(pi/2-DF$dir)*DF$Lmn*f*0.5
  #Compression left point X
  DF$Clef_x<-DF$x - f_com_sin
  #Compression right point X
  DF$Crig_x<-DF$x + f_com_sin
  #Compression left point Y
  DF[idxpos,'Clef_y']<-DF[idxpos,'y'] - f_com_cos[idxpos]
  DF[idxneg,'Clef_y']<-DF[idxneg,'y'] + f_com_cos[idxneg]
  #Compression left point Y
  DF[idxpos,'Crig_y']<-DF[idxpos,'y'] + f_com_cos[idxpos]
  DF[idxneg,'Crig_y']<-DF[idxneg,'y'] - f_com_cos[idxneg]
  
  DF<-DF[complete.cases(DF),]
  
  return(DF)
}



################################################################################
# CONVERT DEGREES TO RADIANS OR RADIANS TO DEGREES
#			  
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}




################################################################################
# SHOW MEMORY USAGE OF SINGLE OBJECTS IN CURRENT WORKSPACE
#
# based on postings by Petr Pikal and David Hinds to the r-help list in 2004
# modified by: Dirk Eddelbuettel (http://stackoverflow.com/questions/1358003/
# tricks-to-manage-the-available-memory-in-an-r-session)
.ls. <- function (pos = 1, pattern, order.by = "Size",
                  decreasing=TRUE, head = TRUE, n = 10)
{
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size) / 10^6 # megabytes
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}						  


##########################################################################
# GENERATE OUTLINES AROUND CORRECTED LOOSE CLIFF PIXELS
#  & GIVE INITIAL ID FROM POINTS TO POLYGONS
#
# pts:             SpatialPoints of geom. corrected cliff pixels 
# resol:           Resolution of DEM
# finalBuffer_sz:  Final (neg.) buffer in which defines new cliff outlines

GenerateOutlines<-function(pts,resol,finalBuffer_sz){
  
  ## for testing only:
  # pts<-pts2_spdf
  # resol<-resol
  
  ifelse(resol > 1,
         buffer<-resol,
         buffer<-1
  )
  
  # create outside buffers around each point(width=res)
  bufOut<-gBuffer(pts,width=buffer,byid=TRUE)
  
  # unify circles where they intersect
  bufOut2<-gUnaryUnion(bufOut,id=bufOut$ID) 
  
  # create additional buffer around unified polygon(s)
  #  to reduce holes etc.
  bufOut3<-gBuffer(bufOut2,width=resol/2)  
  
  ## reduce polygon shape with neg. (inside) buffer (width=resol)
  ##  (then the final polygon covers +/- the area of all new cliff pixels) 
  new_cliffs_p<-gBuffer(bufOut3,width=-finalBuffer_sz,byid=TRUE)
  
  ## check:
  # e<-bboxM(cliffs_p,bufOut3)
  # e<-as(extent(e),'SpatialPolygons') 
  # projection(e)<-projec  
  # plot(e,col='white',border='NA')
  # plot(cliffs_p,col='firebrick3',border=NA,add=TRUE)
  # plot(bufOut,border='gray50',lwd=2,add=TRUE)
  # plot(bufOut2,border='orange',lwd=5,add=TRUE)
  # plot(bufOut3,border='dodgerblue',lwd=4,add=TRUE)
  # plot(new_cliffs_p,border='forestgreen',lwd=4,add=TRUE)
  # plot(pts,lwd=3,col='black',add=TRUE)	
  
  # problems with 'gBuffer' evt.?
  # https://gis.stackexchange.com/questions/177638/gbufferrgeos-causing-r-to-crash
  
  plot(new_cliffs_p)  
  ## check:	
  #ids<-sapply(slot(new_cliffs_p,'polygons'), function(x) slot(x,'ID'))
  #plot(new_cliffs_p)
  #text(coordinates(new_cliffs_p),labels=ids,cex=0.7)
  
  # still remaining holes inside the new cliff polygons can be removed  
  p<-slot(new_cliffs_p,'polygons')
  holes<-lapply(p, function(x) sapply(slot(x,'Polygons'),slot,'hole'))
  filled<-lapply(1:length(p),
                 function(i) slot(p[[i]],'Polygons')[!holes[[i]]])
  ids<-row.names(new_cliffs_p)
  new_cliffs_fill_p<-SpatialPolygons(lapply(1:length(filled), 
                                            function(i) Polygons(filled[[i]],ID=ids[i])),
                                     proj4string=CRS(proj4string(new_cliffs_p)))
  
  ## check:
  #plot(new_cliffs_fill_p,col='red')
  #plot(cliffs_p,lwd=1,add=TRUE)
  #plot(pts2,col='white',add=TRUE)
  rm(new_cliffs_p)
  return(new_cliffs_fill_p)  
}


###########################################################
# DETECT & REMOVE TRANSITION ZONE BETWEEN OLD AND NEW CLIFF
#
# xyz_m:          Old matrix with x/y/elevation values of initial cliff
# new_xyz_m:      New matrix with x/y/elevation values of melted cliff
# DEM:            DEM-raster
# LkAffected_id:  Vector of 1/0 binary code for lake-affected/not-affected
#                  old cliff pixels

RemoveTransZone<-function(xyz_m,new_xyz_m,dem,LkAffected_id){
  
  
  ## only for testing:
  # xyz_m<-extd_xyz_m
  # new_xyz_m<-new_xyz_m
  # dem<-gl_geom_st$Elevation
  # LkAffected_id<-LkAffected_id
  
  
  # create tracks (lines) leading from old to new cliff pixels
  tracks<-list()
  for (irow in 1:nrow(xyz_m)) {
    tracks[[irow]]<-Lines(Line(rbind(as.numeric(xyz_m[irow,c('x','y')]),
                                     as.numeric(new_xyz_m[irow,c('x','y')]))),
                          ID=as.character(irow))
  }
  tracks<-SpatialLines(tracks)
  projection(tracks)<-projection(dem)      
  
  ## only for testing:
  # plot(dem)
  # Arrows(xyz_m[,'x'],xyz_m[,'y'],new_xyz_m[,'x'],new_xyz_m[,'y'],code=2,
  # arr.length=0.1,arr.adj=1,arr.type='triangle')
  ## e<-drawExtent()
  # plot(crop(dem,e))      #plot(dem)
  # Arrows(xyz_m[,'x'],xyz_m[,'y'],new_xyz_m[,'x'],new_xyz_m[,'y'],code=2,
  # arr.length=0.1,arr.adj=1,arr.type='triangle')
  
  # extract raster cells by direct line from old to new cliff pixel
  #  (output: cellnumber & elevation)
  linecells_list<-extract(dem,tracks,cellnumbers=TRUE) # !! SLOW !!
  
  # number of raster cells selected by each line from old to new cliff pixel
  length_v<-sapply(linecells_list,nrow)
  
  # extract elevation values of extended old cliff cells from DEM    
  elevation_extd<-extract(dem,xyz_m)
  
  
  #    # decrease in elevation [m] to be applied to raster cells between 
  #    #  old & new cliff pixels
  #    elevationseq_list<-list()
  #    for(i in 1:nrow(xyz_m)){
  #        el_seq<-seq(xyz_m[i,'elevation'],new_xyz_m[i,'elevation'],
  #                    length=length_v[i])
  #        elevationseq_list[[i]]<-el_seq
  #        }
  
  #index with all old coordinates within lake-affected area
  idx<-which(LkAffected_id == 1)
  
  # decrease in elevation [m] to be applied to raster cells between 
  #  old & new cliff pixels    
  elevationseq_list<-list()
  
  
  for(i in 1:nrow(xyz_m)){
    el_seq<-seq(elevation_extd[i],new_xyz_m[i,'elevation'],
                length=length_v[i])   
    # if old cliff pixel is lake-affected, keep elevation constant
    if(i %in% idx){el_seq[1:length(el_seq)]<-el_seq[1]}
    elevationseq_list[[i]]<-el_seq
  }
  
  # extract all values from DEM-raster
  DEM_vals<-extract(dem,1:ncell(dem))
  origDEM_vals<-DEM_vals   
  DEM_vals<-cbind(replicate(nrow(xyz_m),DEM_vals)) # for each point
  #    DistEle_ls<-list()
  
  # apply new elevation values to raster cells under each line via index
  for(j in 1:nrow(xyz_m)){
    ele<-elevationseq_list[[j]] 
    idx<-linecells_list[[j]][,'cell']
    linecoord_v<-xyFromCell(dem,idx)          
    distToOrigin<-pointDistance(xyz_m[j,1:2],linecoord_v,lonlat=FALSE)
    dfr<-data.frame(cbind(idx,distToOrigin,ele))
    names(dfr)<-c('cell','d','ele')
    # distances from original cells to each cell extracted by line, 
    #  ordered from the closest to the farthest 
    #  (attached: cell number)
    ord_dfr<-dfr[order(dfr$d),]
    #        DistEle<-ord_dfr[,-1]
    #        DistEle$Pt<-'Pt'&rep(j,nrow(DistEle))
    #        DistEle_ls[[j]]<-DistEle           
    #########
    # IMPLEMENT BACKWASTING OVER RIDGE/SADDLE !!
    #
    #
    #########  
    DEM_vals[ord_dfr$cell,j]<-ord_dfr$ele
  }
  
  
  ## testing only:
  #library(plyr)
  #DistEle.long<-melt(ldply(DistEle_ls,data.frame),'Pt')
  
  # set to NA if cell is not changed (not affected by backwasting track)
  # (otherwise apply(...,mean) not meaningfull)
  DEM_vals<-apply(DEM_vals[,1:ncol(DEM_vals)],2,function(x) x-origDEM_vals)
  #    DEM_vals<-data.frame(DEM_vals)
  DEM_vals[DEM_vals == 0]<-NA    
  
  # take always smallest elevation value per affected raster cell
  #    DEM_vals<-apply(DEM_vals,1,min,na.rm=TRUE)
  DEM_vals<-apply(DEM_vals,1,mean,na.rm=TRUE)
  # give DEM back with new elevation values for affected raster cells
  idx<-which(is.na(DEM_vals))                       #find cells not affected
  DEM_vals[-idx]<-origDEM_vals[-idx]+DEM_vals[-idx] #transition cells
  DEM_vals[idx]<-origDEM_vals[idx]                  #not affected cells
  new_DEM<-setValues(dem,DEM_vals)
  
  ## only for testing:
  #    plot(crop(dem,e))
  #    plot(crop(new_DEM,e))
  #    plot(new_DEM-dem)    
  #    plot(crop(new_DEM-dem,e))
  #    Arrows(xyz_m[,'x'],xyz_m[,'y'],new_xyz_m[,'x'],new_xyz_m[,'y'],code=2,
  #           arr.length=0.1,arr.adj=1,arr.type='triangle')
  #    summary(new_DEM-dem)
  
  return(new_DEM)
}


################################################################################
# MODE
#
# v:          Vector with values

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


################################################################################
# R SQUARED
#
# x:          Vector with values
# y:          Vector with values

rsq <- function (x,y) cor(x,y,use='complete.obs') ^ 2




################################################################################
# CIRCULAR MEAN (OF ASPECTS/DIRECTIONS)
#
# aspects_v:          Vector with aspects [?] per pixel to average

# https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf
# To interpret these quantities it may be useful to imagine that each angle 
# represents a vector of length one in the direction of the angle. 
# Suppose these individual vectors are arranged so that the beginning of the 
# first vector is at the origin, the beginning of the second vector is at the
# end of the first, the beginning of the third vector is at the end of the 
# second, and so on. We can then imagine a single vector a that will stretch 
# from the origin to the end of the last observation. 

CircMean<-function(aspects_v){
  
  rads<-aspects_v*pi/180    # turn values from deg. into rad.
  avsin<-mean(sin(rads),na.rm=T)
  avcos<-mean(cos(rads),na.rm=T)
  avaz<-atan2(avsin,avcos)
  if (!is.na(avaz) & avaz<0){avaz<-avaz+(2*pi)}
  meanAsp<-avaz*180/pi
  return(meanAsp)
}



################################################################################
# CIRCULAR VARIANCE (OF ASPECTS/DIRECTIONS)
#
# aspects_v:          Vector with aspects [?] per pixel to average

# https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf

CircVar<-function(aspects_v){
  
  rads<-aspects_v*pi/180    # turn values from deg. into rad.
  ssin<-sum(sin(rads))
  scos<-sum(cos(rads))
  rbar<-sqrt(ssin^2+scos^2)/length(rads)
  varAsp<-1-rbar
  return(varAsp)
  # Note that varAsp varies between 0 and 1 and that 
  # a value of varAsp near 0 implies that there was 
  # little variation in values of the angles.
}



################################################################################
# CIRCULAR STANDARD DEVIATION (OF ASPECTS/DIRECTIONS)
#
# aspects_v:          Vector with aspects [?] per pixel to average

# https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf

CircSD<-function(aspects_v){
  
  rads<-aspects_v*pi/180    # turn values from deg. into rad.
  ssin<-sum(sin(rads))
  scos<-sum(cos(rads))
  rbar<-sqrt(ssin^2+scos^2)/length(rads)
  sdAsp<-sqrt(-2*log(rbar))
  sdAsp<-sdAsp*180/pi
  return(sdAsp)
}



################################################################################
# WEIGHTED CIRCULAR MEAN (OF ASPECTS/DIRECTIONS)
#
# aspects_v:          Vector with aspects  [?] per pixel to average
# weight_v:           Vector with weights (e.g. inclined area) per pixel 
#                      to average

CircMeanW<-function(aspects_v,weight_v){
  
  ## FOR TESTING:
  #aspects_v<-c(90,270,270,180)
  #weight_v<-c(1.5,1,1,1.5)
  
  rads<-aspects_v*pi/180    # turn values from deg. into rad.
  avsin<-sum(sin(rads)*weight_v)/sum(weight_v)
  avcos<-sum(cos(rads)*weight_v)/sum(weight_v)
  avaz<-atan2(avsin,avcos)
  if (avaz<0){avaz<-avaz+(2*pi)}
  meanAsp<-avaz*180/pi
  return(meanAsp)
}



################################################################################
# AGGREGATION OF VECTORS (E.G. ACCORDING TO DATES) TO WEIGHTED MEANS
#
# DATA:      Vector containing data to aggregate by weighted mean 
# ID:        Vector containing identifier for aggregation (e.g. "Date")
# W:         Vector containing weights per element 
# CircMean   Define (TRUE or FALSE) if vectorial mean has to be calculated

AggregateWmean<-function(DATA,ID,W,CircMean){
  
  ## FOR TESTING:
  # DATA<-df$Aspect
  # ID<-df$Date
  # W<-df$w
  # CircMean<-FALSE
  
  df<-as.data.frame(cbind(DATA,ID,W))
  if(CircMean==FALSE){               
    wMean<-as.vector(by(df[c('DATA','W')],
                        list(df$ID),
                        function(x){do.call(weighted.mean,unname(x))}
    )
    )
  }
  
  #TODO: CHECK VECTORIAL MEAN FUNCTION AGAIN
  if(CircMean==TRUE){               
    wMean<-as.vector(by(df[c('DATA','W')],
                        list(df$ID),
                        function(x){do.call(CircMeanW,unname(x))}
    )
    )
  }				
  return(wMean)
}



################################################################################
# ADDING AN (JPG, PNG) IMAGE TO A PLOT
#
# logo:   image (read to R-session using "readJPEG()")
# x, y:   fraction of the the usr space where image should be placed
# size:   relative size (respects aspect ratio)
#
# https://gist.github.com/scrogster/7fc5b7597b63585a00b6

add_img<-function(logo, x, y, size){
  
  dims<-dim(logo)[1:2]#number of x-y pixels for the image (aspect ratio)
  AR<-dims[1]/dims[2]
  par(usr=c(0, 1, 0, 1))
  rasterImage(logo,
              x-(size/2),
              y-(AR*size/2),
              x+(size/2),
              y+(AR*size/2),
              interpolate=TRUE
  )
}



################################################################################
# (MATRIX-) ROTATION OF A RASTER AROUND A CENTRE POINT
#
# theta:    Rotation angle [?] (pos.=clockwise; neg.=anti-clockwise)
# r:        Raster to rotate [RasterLayer]
# cpt:      Rotation centre point [SpatialPoints] within raster extent
# resol:    Resolution of original raster [m]

RasterRotation<-function(theta_deg,cpt,r,resol){
  
  plot(r)                     
  
  # create xyz-dataframe from raster
  r_df<-as.data.frame(coordinates(r))
  r_df$z<-as.data.frame(r)[,1]
  
  # conversion of rotation angle to fulfill: 
  #  pos.=clockwise; neg.=anti-clockwise
  theta_deg<-theta_deg*-1
  # convert rotation angle into radians
  theta<-theta_deg*pi/180
  
  # calculate new coordinates after rotation
  r_df$Xr<-r_df$x * cos(theta) - r_df$y * sin(theta)
  r_df$Yr<-r_df$x * sin(theta) + r_df$y * cos(theta)
  
  # get centre coordinates
  ccoord<-coordinates(cpt)
  ccoord<-as.data.frame(extract(r,ccoord,cellnumbers=TRUE))
  ccoord<-r_df[ccoord[,1],1:2]
  id<-as.numeric(rownames(ccoord))
  dCx<-r_df[id,'x']-r_df[id,'Xr']
  dCy<-r_df[id,'y']-r_df[id,'Yr']
  
  # shift rotated raster back to rotation point
  r_df$Xr<-r_df$Xr+dCx
  r_df$Yr<-r_df$Yr+dCy
  
  # create raster for rotated xyz-dataframe                
  r_m<-as.matrix(r_df)
  ext<-extent(r_m[,4:5])
  res<-res(r)[1]
  r2<-raster(ext,ncol=(ext[2]-ext[1])/res,
             nrow=(ext[4]-ext[3])/res)
  projection(r2)<-projection(r)
  
  # remove NA-cells (NA-cells occur where glacier mask is crossed)
  idx<-which(is.na(r_df$z))
  if(length(idx) != 0){r_df<-r_df[-idx,]}
  
  # linear interpolation (Akima) of raster values 
  #  as data gaps occured due to rotation (->resolution issue)
  xlen2<-ncol(r2)
  ylen2<-nrow(r2)
  xmn2<-mround(extent(r2)[1]+res/2,resol)
  xmx2<-mround(extent(r2)[2]-res/2,resol)
  ymn2<-mround(extent(r2)[3]+res/2,resol)
  ymx2<-mround(extent(r2)[4]-res/2,resol)
  #                xmn2<-extent(r2)[1]-res/2
  #                xmx2<-extent(r2)[2]-res/2
  #                ymn2<-extent(r2)[3]+res/2
  #                ymx2<-extent(r2)[4]-res/2
  lin_val<-interp(r_df[,'Xr'],r_df[,'Yr'],r_df[,'z'],
                  duplicate='user',dupfun='mean',linear=TRUE,
                  extrap=FALSE, 
                  xo=seq(xmn2,xmx2,length=xlen2),
                  yo=seq(ymn2,ymx2,length=ylen2))
  
  r2<-raster(lin_val)
  projection(r2)<-as.character(projection(r))
  plot(r2)
  
  return(r2)
  rm(r_df,theta,ccoord,id,dCx,dCy,r_df,r_m,ext,res,idx,
     xlen2,ylen2,xmn2,xmx2,ymn2,ymx2,lin_val,r2)
}





################################################################################
# (MATRIX-) ROTATION OF A POLYGON AROUND A CENTRE POINT
#
# theta:    Rotation angle [?] (pos.=clockwise; neg.=anti-clockwise)
# p:        Polygon to rotate [SpatialPolygon]
# cpt:      Rotation centre point [SpatialPoints]

PolygonRotation<-function(theta_deg,cpt,p){
  
  plot(p)     
  
  # create xy-dataframe from polygon
  p_df<-as.data.frame(p@polygons[[1]]@Polygons[[1]]@coords )
  colnames(p_df)<-c('x','y')
  
  # conversion of rotation angle to fulfill: 
  #  pos.=clockwise; neg.=anti-clockwise
  theta_deg<-theta_deg*-1
  # convert rotation angle into radians
  theta<-theta_deg*pi/180
  
  # calculate new coordinates after rotation
  p_df$Xr<-p_df$x * cos(theta) - p_df$y * sin(theta)
  p_df$Yr<-p_df$x * sin(theta) + p_df$y * cos(theta)
  
  # correct for rotation around central point:
  # check
  plot(p)
  plot(cpt,add=TRUE)
  # get centre coordinates
  ccoord<-coordinates(cpt)
  
  ccoordXr<-ccoord[1] * cos(theta) - ccoord[2] * sin(theta)
  ccoordYr<-ccoord[1] * sin(theta) + ccoord[2] * cos(theta)
  
  dCx<-ccoord[1]-ccoordXr
  dCy<-ccoord[2]-ccoordYr
  
  # shift rotated raster back to rotation point
  p_df$Xr<-p_df$Xr+dCx
  p_df$Yr<-p_df$Yr+dCy
  
  # create polygon from rotated vertices
  vertices<-cbind(p_df$Xr,p_df$Yr)
  p_ls<-list()
  p_ls[[1]]<-Polygon(vertices)   
  p2<-Polygons(p_ls,'Poly')
  p2<-SpatialPolygons(list(p2))
  projection(p2)<-as.character(projection(p))
  
  plot(p2,add=TRUE)
  
  return(p2)
  rm(p_df,theta,ccoord,ccoordXr,ccoordYr,dCx,dCy,vertices,p_ls,p2)
}



################################################################################
# (MATRIX-) ROTATION OF A LINE AROUND A CENTRE POINT
#
# theta:    Rotation angle [?] (pos.=clockwise; neg.=anti-clockwise)
# ln:       Line to rotate [SpatialLine]
# cpt:      Rotation centre point [SpatialPoints]

LineRotation<-function(theta_deg,cpt,ln){
  
  plot(ln)     
  
  # create xy-dataframe from polygon
  ln_df<-as.data.frame(coordinates(ln[1,])[[1]][[1]])
  colnames(ln_df)<-c('x','y')
  
  # conversion of rotation angle to fulfill: 
  #  pos.=clockwise; neg.=anti-clockwise
  theta_deg<-theta_deg*-1
  # convert rotation angle into radians
  theta<-theta_deg*pi/180
  
  # calculate new coordinates after rotation
  ln_df$Xr<-ln_df$x * cos(theta) - ln_df$y * sin(theta)
  ln_df$Yr<-ln_df$x * sin(theta) + ln_df$y * cos(theta)
  
  # correct for rotation around central point:
  # check
  plot(ln)
  plot(cpt,add=TRUE)
  # get centre coordinates
  ccoord<-coordinates(cpt)
  
  ccoordXr<-ccoord[1] * cos(theta) - ccoord[2] * sin(theta)
  ccoordYr<-ccoord[1] * sin(theta) + ccoord[2] * cos(theta)
  
  dCx<-ccoord[1]-ccoordXr
  dCy<-ccoord[2]-ccoordYr
  
  # shift rotated raster back to rotation point
  ln_df$Xr<-ln_df$Xr+dCx
  ln_df$Yr<-ln_df$Yr+dCy
  
  # create line from rotated line start/end points
  ln2<-SpatialLines(list(Lines(Line(cbind(ln_df$Xr,ln_df$Yr)),
                               ID='a')))
  projection(ln2)<-as.character(projection(ln2))
  
  plot(ln2,add=TRUE)
  plot(ln2)
  plot(cpt,add=TRUE)
  
  return(ln2)
  rm(ln_df,theta,ccoord,ccoordXr,ccoordYr,dCx,dCy,vertices,ln2)
}



################################################################################
# SUPERIMPOSE DEM:
#  
#  Superimpose a DEM over another DEM for cases merge() or mosaic() do not work 
#   (e.g. when small inaccuracies occur due to raster rotation).
#   Merge is based on coordinates only.
#
#  topDEM:    DEM to superimpose (DEM put on top)
#  baseDEM:   DEM (serving as base) being overlayed by topDEM 


SuperimposeDEM <- function(topDEM,baseDEM){
  
  ele_clip<-as.data.frame(topDEM)
  coord_clip<-coordinates(topDEM)
  xyz_clip<-cbind(coord_clip[,'x'],coord_clip[,'y'],ele_clip)
  colnames(xyz_clip)<-c('x','y','z')
  xyz_clip<-xyz_clip[complete.cases(xyz_clip),]
  
  # extract cellnumbers of all raster pixels in the base DEM 
  #  corresponding to the top DEM cells 
  #  (necessary step due to coordinate inaccuraciees)
  idx<-extract(baseDEM,xyz_clip[,1:2],cellnumbers=TRUE)[,1]
  allcells<-as.vector(baseDEM)
  allcells[idx]<-xyz_clip$z
  newDEM<-setValues(baseDEM,allcells)
  return(newDEM)
}



################################################################################
# MERGING MULTIPLE DATA FILES INTO ONE DATA FRAME
#  http://stackoverflow.com/questions/8091303/simultaneously-merge-
#  multiple-data-frames-in-a-list
#
# df_ls:          list of dataframes, all having one column name in common
# MatchString:    string specifing column name to match

MultiMerge<-function(df_ls,MatchString){
  
  # to avoid (warning message about) duplicated column names:
  #  ->make columns to merge unique by adding a number
  df_ls<-Map(function(x, i) setNames(x,ifelse(names(x) %in% MatchString,
                                              names(x),sprintf('%s.%d',names(x),i))),df_ls,
             seq_along(df_ls)
  )
  
  # merge all dataframes based on single column ("MatchString")        
  Reduce(function(x,y) merge(x,y,MatchString,all=T,suffixes=seq(1:6)),
         df_ls)
}


################################################################################
# CALCULATE DIURNAL CYCLE (+SD) PER CLIFF FROM MATRIX ->TS x CLIFF PIXELS
#
# cliff_IDs:   Vector with cliff ID per pixel
# ts_v:        Vector with timesteps (numeric or POSIX-format) 
# Flux_m:      Matrix with rows as timesteps (equal to "ts_v") 
#               and columns for each cliff pixel (of different cliffs)

DC_PC<-function(cliff_IDs,ts_v,Flux_m){
  
  unique_ids<-unique(cliff_IDs)
  date_h<-as.POSIXlt(ts_v,origin='1970-01-01')$hour
  
  Flux_dc_ls<-list()
  Flux_sd_ls<-list()
  for(cl in 1:length(unique_ids)){
    # create index to select all pixels of a single cliff
    idx<-which(cliff_IDs == unique_ids[cl])
    
    # aggregate flux values to diurnal cycle per cliff pixel            
    Flux_dc<-aggregate(Flux_m[,idx],list(date_h),mean)
    
    # derive one diurnal cycle for entire cliff
    Flux_dc_ls[[cl]]<-apply(Flux_dc,1,mean)
    
    # calculate standard deviation between all pixels of a cliff 
    #  per hour of diurnal cycle
    Flux_sd_ls[[cl]]<-apply(Flux_dc,1,sd)
  }
  output_ls<-list(Flux_dc_ls,Flux_sd_ls)
  names(output_ls)<-c('DC','SD')
  return(output_ls)
}



################################################################################
DC<-function(dataframe){
  
  # CALCULATE DIURNAL CYCLE (+SD) PER DATA-VECTOR ->TS x VARIABLES
  #
  # dataframe:   Dataframe with timesteps (numeric or POSIX-format) 
  #               & values for each timestep (each column one variable)
  #               ->only one column in POSIXct-format allowed (date-column)
  #
  # e.g.:
  #  Date                  L_s   L_d   H      Q_m   melt
  #  2013-05-19 00:00:00   220.3 127.5 6.793  48.44 0.0006
  #  2013-05-19 01:00:00   216.2 126.2 3.954  40.19 0.0005
  #  2013-05-19 02:00:00   210.8 122.5 2.172  29.27 0.0004
  #  2013-05-19 03:00:00   208.3 119.5 2.453  24.14 0.0003
  #  2013-05-19 04:00:00   210.6 120.2 3.912  28.48 0.0004
  #  2013-05-19 05:00:00   212.3 121.6 4.608  32.30 0.0004
  idx<-which(sapply(dataframe,is.POSIXct))
  colnames(dataframe)[idx]<-'Date' 
  
  # extract only hours per timestep
  date_h<-as.POSIXlt(dataframe$Date,origin='1970-01-01')$hour
  
  # remove 'Date'-column
  dataframe<-dataframe[!colnames(dataframe) %in% 'Date']
  
  # aggregate flux values to diurnal cycle (and compute sd)
  Flux_dc_ls<-list()
  Flux_sd_ls<-list()
  for(v in 1:ncol(dataframe)){          
    Flux_dc_ls[[v]]<-aggregate(dataframe[[v]],list(date_h),mean,na.rm=TRUE)[,2]
    Flux_sd_ls[[v]]<-aggregate(dataframe[[v]],list(date_h),sd,na.rm=TRUE)[,2]
  }
  
  # diurnal cycle                                       
  DC_df<-data.frame(1:24,do.call(cbind,Flux_dc_ls))
  colnames(DC_df)<-c('hours',colnames(dataframe)) 
  
  # standard deviation
  SD_df<-data.frame(1:24,do.call(cbind,Flux_sd_ls))
  colnames(SD_df)<-c('hours',colnames(dataframe)) 
  
  output_ls<-list(DC_df,SD_df)
  names(output_ls)<-c('DC','SD')
  return(output_ls)
}







################################################################################
# BOUNDING BOX AROUND MULTIPLE RASTERS OR POLYGONS
#
# ... :          Rasters or polygons, connected with c() or list
#
# Output "e":    extent of boundary box around input objects

bboxM <- function(data){
  
  ifelse(is.list(data),
         spobj_ls<-data,
         spobj_ls<-list(data)
  )
  xmn_v<-vector()
  xmx_v<-vector()
  ymn_v<-vector()
  ymx_v<-vector()
  for(i in 1:length(spobj_ls)){
    xmn_v[i]<-extent(spobj_ls[[i]])[1]
    xmx_v[i]<-extent(spobj_ls[[i]])[2]
    ymn_v[i]<-extent(spobj_ls[[i]])[3]
    ymx_v[i]<-extent(spobj_ls[[i]])[4]
  }
  e<-extent(c(min(xmn_v),max(xmx_v),min(ymn_v),max(ymx_v)))  
  return(e)  
}




################################################################################
# ROUNDING FUNCTIONS
#
# function: ROUND x to neareast of specified number y
# -->mround(x,y)
mround <- function(x,base){
  base*round(x/base)
}
# function: ROUND DOWN x to neareast of specified number y
# -->mrounddown(x,y)
mroundd <- function(x,base){
  base*floor(x/base)
}

# function: ROUND UP x to neareast of specified number y
# -->mrounddown(x,y)
mroundu <- function(x,base){
  base*ceiling(x/base)
}




################################################################################
# LINEAR RASTER INTERPOLATION
#
# r:   raster

InterpRaster <- function(r) {
  
  ## create xyz-dataframe from raster
  r_df<-as.data.frame(coordinates(r))
  r_df$z<-as.data.frame(r)[,1]
  
  # remove NA-cells
  idx<-which(is.na(r_df$z))
  if(length(idx) != 0){r_df<-r_df[-idx,]}
  
  ## linear interpolation (Akima) of raster values 
  # source(path_func)
  resol<-res(r)[1]
  xlen<-ncol(r)
  ylen<-nrow(r)
  xmn<-mround(extent(r)[1]+resol/2,resol)
  xmx<-mround(extent(r)[2]-resol/2,resol)
  ymn<-mround(extent(r)[3]+resol/2,resol)
  ymx<-mround(extent(r)[4]-resol/2,resol)
  lin_val<-interp(r_df[,'x'],r_df[,'y'],r_df[,'z'],
                  duplicate='user',dupfun='mean',linear=TRUE,
                  extrap=FALSE, 
                  xo=seq(xmn,xmx,length=xlen),
                  yo=seq(ymn,ymx,length=ylen))
  
  r2<-raster(lin_val)
  projection(r2)<-as.character(projection(r))
  return(r2)
}


################################################################################
# MAXIMUM  VALUE PER COLUMN 
#
# x:   data frame

colMax <- function(x) sapply(x, max, na.rm = TRUE)





################################################################################
# LOAD ALL KML-LAYERS FROM .KML-FILE 
#
# kmlfile:   e.g. 'C:/Documents/fieldsites.KML'

allKmlLayers <- function(kmlfile){
  lyr <- ogrListLayers(kmlfile)
  mykml <- list()
  for (i in 1:length(lyr)) {
    mykml[i] <- readOGR(kmlfile,lyr[i])
  }
  names(mykml) <- lyr
  return(mykml)
}




################################################################################
# FOCAL FOR RASTERBRICK OR RASTERSTACK 
# x:   rasterbrick or rasterstack

multiFocal <- function(x, w=matrix(1, nr=3, nc=3), ...) { 
  
  if(is.character(x)) { 
    x <- brick(x) 
  } 
  # The function to be applied to each individual layer 
  fun <- function(ind, x, w, ...){ 
    focal(x[[ind]], w=w, ...) 
  } 
  
  n <- seq(nlayers(x)) 
  list <- lapply(X=n, FUN=fun, x=x, w=w, ...) 
  
  out <- stack(list) 
  return(out) 
} 



################################################################################
# REMOVE EMPTY LIST FROM LIST
# 
# delete null/empty entries in a list
#  https://stat.ethz.ch/pipermail/r-help/2006-August/111896.html

delete.NULLs<-function(x.list){   
  x.list[unlist(lapply(x.list, length) != 0)]
}



################################################################################
# UNLOAD A PACKAGE
#
# It is possible to have multiple versions of a package loaded at once
#  (for example, if you have a development version and a stable version
#  in different libraries).
#  To detach, guaranteed that all copies are detached, use this function:
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}


###########################################################
# DIVERGE0
#
# Plotting a raster with the color ramp diverging around zero
#  http://stackoverflow.com/questions/33750235/plotting-a-raster-
#  with-the-color-ramp-diverging-around-zero

diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character 
  #       vector of colour names to interpolate, or a colorRampPalette.
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
  p
}


################################################################################
# RASTER TO POLYGON
#
# https://gist.github.com/johnbaums/26e8091f082f2b3dd279
#
# Using gdal_polygonize.py
#
# Convert raster data to a ESRI polygon shapefile and (optionally) 
#  a SpatialPolygonsDataFrame
#
# Note that the following takes an argument outshape, which specifies the 
# path where the converted polygon shapefile will reside. Default is NULL,
# which saves the shapefile to a temporary file. Either way, the function
# returns a SpatialPolygonsDataFrame representation of the shapefile,
# unless readpoly is FALSE, in which case the function returns NULL.

# Define the function
polygonizer <- function(x, outshape=NULL, pypath=NULL, readpoly=TRUE, 
                        fillholes=FALSE, aggregate=FALSE, 
                        quietish=TRUE) {
  # x: an R Raster layer, or the file path to a raster file recognised by GDAL 
  # outshape: the path to the output shapefile (if NULL, a temporary file will 
  #           be created) 
  # pypath: the path to gdal_polygonize.py or OSGeo4W.bat (if NULL, the function 
  #         will attempt to determine the location)
  # readpoly: should the polygon shapefile be read back into R, and returned by
  #           this function? (logical) 
  # fillholes: should holes be deleted (i.e., their area added to the containing
  #            polygon)
  # aggregate: should polygons be aggregated by their associated raster value?
  # quietish: should (some) messages be suppressed? (logical)
  if (isTRUE(readpoly) || isTRUE(fillholes)) require(rgdal)
  if (is.null(pypath)) {
    cmd <- Sys.which('OSGeo4W.bat')
    pypath <- 'gdal_polygonize'
    if(cmd=='') {
      cmd <- 'python'
      pypath <- Sys.which('gdal_polygonize.py')
      if (!file.exists(pypath)) 
        stop("Could not find gdal_polygonize.py or OSGeo4W on your system.") 
    }
  }
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists)) 
      stop(sprintf('File already exists: %s', 
                   toString(paste(outshape, c('shp', 'shx', 'dbf'), 
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.tif')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  
  system2(cmd, args=(
    sprintf('"%s" "%s" %s -f "ESRI Shapefile" "%s.shp"', 
            pypath, rastpath, ifelse(quietish, '-q ', ''), outshape)))
  
  if(isTRUE(aggregate)||isTRUE(readpoly)||isTRUE(fillholes)) {
    shp <- readOGR(dirname(outshape), layer=basename(outshape), 
                   verbose=!quietish)    
  } else return(NULL)
  
  if (isTRUE(fillholes)) {
    poly_noholes <- lapply(shp@polygons, function(x) {
      Filter(function(p) p@ringDir==1, x@Polygons)[[1]]
    })
    pp <- SpatialPolygons(mapply(function(x, id) {
      list(Polygons(list(x), ID=id))
    }, poly_noholes, row.names(shp)), proj4string=CRS(proj4string(shp)))
    shp <- SpatialPolygonsDataFrame(pp, shp@data)
    if(isTRUE(aggregate)) shp <- aggregate(shp, names(shp))
    writeOGR(shp, dirname(outshape), basename(outshape), 
             'ESRI Shapefile', overwrite=TRUE)
  }
  if(isTRUE(aggregate) & !isTRUE(fillholes)) {
    shp <- aggregate(shp, names(shp))
    writeOGR(shp, dirname(outshape), basename(outshape), 
             'ESRI Shapefile', overwrite=TRUE)
  }
  ifelse(isTRUE(readpoly), return(shp), return(NULL))
}



################################################################################
# PLOT SETTINGS:
# Custom themes for ggplot2
#
# see example:
#  http://stackoverflow.com/questions/6736378/how-do-i-change-the-background-
#  color-of-a-plot-made-with-ggplot2

# raster-plots (single, black background)
theme_rasterSingleBlack<-function(base_size=15,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      panel.background=element_rect(fill='black'),
      panel.border=element_blank(),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(colour='grey60')
    )
}


################################################################################
# raster-plots (single, light grey background)
theme_rasterSingleBright<-function(base_size=15,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      panel.background=element_rect(fill='grey90'),
      panel.border=element_blank(),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(colour='grey60')
    )
}


################################################################################
# raster-plots (single, white background)
theme_rasterSingleWhite<-function(base_size=9,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      legend.text=element_text(size=9),
      legend.title=element_text(size=9),
      axis.text=element_text(colour='black'),
      axis.ticks=element_blank(),
      panel.border=element_rect(colour='white',fill=NA,size=0.5),
      panel.background=element_rect(fill='grey90'),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(colour='white',size=0.5)
    )
}


################################################################################
# raster-plots (single, white background, no axes)
theme_rasterSingleWhite_noAxes<-function(base_size=9,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      axis.ticks=element_blank(),
      legend.text=element_text(size=9),
      legend.title=element_text(size=9),
      panel.border=element_rect(colour='white',fill=NA,size=0.5),
      panel.background=element_rect(fill='grey90'),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(colour='white',size=0.5)
    )
}


################################################################################
# raster-plots (multiple)
theme_rasterMultiple<-function(base_size=15,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      plot.background = element_blank(),
      panel.background=element_rect(fill='white'),
      panel.grid.minor = element_blank(),
      panel.grid.major=element_line(colour='grey60'),
      panel.border = element_blank(),
      strip.background=element_blank(),
      axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank()
    )
}


################################################################################
# line-plots (legend on top)
theme_line<-function(base_size=15,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      legend.title = element_blank(),
      legend.position = 'top',
      legend.direction = 'horizontal',
      axis.ticks = element_blank(),
      axis.text = element_text(colour='black')
    )
}


################################################################################
# line-plots
theme_line2<-function(base_size=15,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      legend.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(colour='black')
    )
}


################################################################################
# line-plots (legend right)
theme_line3<-function(base_size=15,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      legend.title = element_blank(),
      legend.position = 'top',
      legend.direction = 'horizontal',
      axis.ticks = element_blank(),
      axis.text = element_text(colour='black')
    )
}



################################################################################
# line-plots (no legend)
theme_line4<-function(base_size=15,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      legend.title = element_blank(),
      legend.position = 'none',
      axis.ticks = element_blank(),
      axis.text = element_text(colour='black')
    )
}


################################################################################
# bar-plots
theme_bar<-function(base_size=20,base_family='Arial'){
  theme_gray(base_size=base_size,base_family=base_family) %+replace%
    theme(
      axis.text = element_text(colour='black'),
      axis.ticks = element_line(colour='grey60',size=0.25),
      panel.background = element_blank(),
      panel.grid.major=element_line(colour='grey60',size=0.25),
      panel.grid.minor=element_line(colour='grey60',size=0.25),
      panel.border = element_blank(),
      legend.text = element_text(size=20)
    )
}



################################################################################
#extract legend: share a legend between two ggplot2 graphs
#https://github.com/hadley/ggplot2/wiki/
# Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp<-ggplot_gtable(ggplot_build(a.gplot))
  leg<-which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
  legend<-tmp$grobs[[leg]]
  return(legend)}



################################################################################
# Multiple plot function
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# e.g. multiplot(p1,p2,cols=1)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



################################################################################
#extract RColorBrewer palette colors
# pal:     palette as character string from RColorBrewer palettes

extrBrewerCol <- function(pal) brewer.pal(brewer.pal.info[pal,'maxcolors'],pal)



################################################################################
#add transparency information to colour string 
# col:     vector of colour coded strings (e.g. in the form "#FF7F00", "...",)
# alpha:   transparency

col2alpha <- function(col, alpha) {
  alphacolor_v<-vector()
  for(i in 1:length(col)){
    col_rgb <- col2rgb(col[i])/255
    alphacolor_v[i]<-rgb(col_rgb[1], col_rgb[2], col_rgb[3], alpha = alpha)
  }
  alphacolor_v
}



################################################################################
#pairwise difference operator 
# https://grokbase.com/t/r/r-help/047ytccbrr/r-pairwise-difference-operator
# m:     matrix
pairwiseDiff <- function(m){
  npairs <- choose( ncol(m), 2 )
  results <- matrix( NA, nc=npairs, nr=nrow(m) )
  cnames <- rep(NA, npairs)
  if(is.null(colnames(m))) colnames(m) <- paste("col", 1:ncol(m), sep="")
  
  k <- 1
  for(i in 1:ncol(m)){
    for(j in 1:ncol(m)){
      if(j <= i) next;
      results[ ,k] <- m[ ,i] - m[ ,j]
      cnames[k] <- paste(colnames(m)[ c(i, j) ], collapse=".vs.")
      k <- k + 1
    }
  }
  
  colnames(results) <- cnames
  rownames(results) <- rownames(m)
  return(results)
}




##############################################################################
# Filename: matlab_time.R
# Convert between MATLAB datenum values and R POSIXt time values.
# 
# Author: Luke Miller   Feb 20, 2011
# https://lukemiller.org/index.php/2011/02/converting-matlab-and-r-date-and-time-values/
#############################################################################

#Convert a numeric  MATLAB datenum (days since 0000-1-1 00:00) to seconds in 
#the Unix epoch (seconds since 1970-1-1 00:00). Specify a time zone if the 
#input datenum is anything other than the GMT/UTC time zone. 
matlab2POS = function(x, timez = "UTC") {
  days = x - 719529 	# 719529 = days from 1-1-0000 to 1-1-1970
  secs = days * 86400 # 86400 seconds in a day
  # This next string of functions is a complete disaster, but it works.
  # It tries to outsmart R by converting the secs value to a POSIXct value
  # in the UTC time zone, then converts that to a time/date string that 
  # should lose the time zone, and then it performs a second as.POSIXct()
  # conversion on the time/date string to get a POSIXct value in the user's 
  # specified timezone. Time zones are a goddamned nightmare.
  
  # return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
  #                                       tz = 'UTC'), format = '%Y-%m-%d %H:%M', 
  #                            tz = 'UTC', usetz = FALSE), tz = timez))
  return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
                                        tz = 'UTC'), 
                             format = '%Y-%m-%d %H:%M:%OS', 
                             tz = 'UTC', usetz = FALSE), tz = timez))
}


###############################################################################
#Convert POSIXct, POSIXlt or 'seconds since 1970-1-1' to MATLAB datenum value.
#The conversion drops any time zone associated with the POSIXt value. It is the
#user's responsibility to keep track of time zones in MATLAB datenums.
#The constant 719529 in the function is the days from 0000-1-1 to 1970-1-1.
POSIXt2matlab = function(x) {
  if (class(x)[1] == "POSIXlt"){
    days = as.numeric(as.Date(x)) #extract days since 1970-1-1
    frac.day = (((x$hour)*3600) + ((x$min)*60) + x$sec)/86400
    datenum = 719529 + days + frac.day 
    datenum = 719529 + days + frac.day		
  } else if (class(x)[1] == "POSIXct"){
    x = as.POSIXlt(x) #convert to POSIXlt class
    days = as.numeric(as.Date(x)) #extract days since 1970-1-1
    frac.day = (((x$hour)*3600) + ((x$min)*60) + x$sec)/86400
    datenum = 719529 + days + frac.day
  } else if (class(x)[1] == "numeric"){
    days = x / 86400 #convert seconds to days
    datenum = days + 719529 
  } else {
    stop("Input cannot be coerced to POSIXlt or numeric value")
  }
  return(datenum)
}
#The output is a numeric vector of 'days since 0000-1-1 00:00'. 



###############################################################################
#Convert POSIXct or POSIXlt objects to MATLAB datenum, in UTC time zone. 
#All time stamps with non-GMT/UTC time zones will be first converted to the 
#GMT/UTC time zone, then converted to MATLAB datenum value. 
POSIXt2matlabUTC = function(x) {
  if (class(x)[1] == "POSIXct") {
    x = as.POSIXlt(x, tz = "UTC") #convert to UTC time zone
    days = as.numeric(x) / 86400 #convert to days
    datenum = days + 719529 #convert to MATLAB datenum
  } else if (class(x)[1] == "POSIXlt") {
    x = as.POSIXlt(x, tz = "UTC") #convert to UTC time zone
    days = as.numeric(x) / 86400  #convert to days
    datenum = days + 719529 #convert to MATLAB datenum
  } else {stop("POSIXct or POSIXlt object required for input")}
  return(datenum)
}
#The output is a numeric vector of 'days since 0000-1-1 00:00', adjusted to UTC



###############################################################################
#DO ROW REPEAT & COL REPEAT
#
# x:  vector to be repeated 
# n:  number of replication
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}






###############################################################################
# PARTITION LIQUID-SOLID PRECIPITATION (DING ET AL., 2014)
#
# Ding et al., 2014 JH  (http://dx.doi.org/10.1016/j.jhydrol.2014.03.038)
#
# Ta:     vector of air temperature [°C]
# rH:     vector of air relative humidity [%]
# Pres:   vector of atmospheric pressure [Pa]

PPARTITION_DING2014 = function(Ta,rH,Pres) {
  
  ## Constants
  g<-9.81 # [m/s^2] gravity acceleration
  P_Ref<-1013.25 # [Pa] reference pressure
  Rd<-287.05 # [J/kg K] dry air gas constant
  
  esat<-611*exp(17.27*Ta/(237.3+Ta)) # [Pa] Vapor pressure saturation
  ea<-esat*rH/100                        #Vapour pressure (Pa) 
  RH<-ea/esat # Relative Humidity
  Laten<-1000*(2501.3 - 2.361*(Ta)) # Latent heat vaporization/condensaition [J/kg]
  cp<-1005 + ((Ta +23.15)^2)/3364 # specific heat air  [J/kg K]
  gam<-cp*100*Pres/(0.622*Laten) # [Pa/C] psycrometric constant
  del<-(4098*esat)/((237.3+Ta)^2) # Pa/C
  Twb<-AWS$Ta - ( esat - ea )/( gam + del) # [C] Wet bulb temperature
  
  ## Reference elevation
  Zref<- -((Ta+15)/2+273.15)*(Rd/g)*log(Pres/P_Ref) # [m]
  Zref<-Zref/1000 # elevation, [m] => [km]
  
  ## Thresholds for discrimination precipitation
  dT<-0.215 - 0.099*RH + 1.018*RH^2
  dS<-2.374 - 1.634*RH
  T0<- -5.87 - 0.1042*Zref + 0.0885*Zref^2 + 16.06*RH - 9.614*RH^2
  Tmin<-T0
  Tmax<-T0
  
  idx<-dT/dS > log(2)
  Tmin[idx]<-T0[idx] - 
    dS[idx]*log(exp(dT[idx]/dS[idx]) - 
                  2*exp(-dT[idx]/dS[idx]))
  Tmax[idx]<-2*T0[idx] - Tmin[idx]
  rm(idx)
  
  ## Calculation of solid fraction of precipitation
  solid_fraction<-1/(1 + exp((Twb - T0)/dS)) # sleet
  
  Pr_sno<-vector()
  Pr_liq<-vector()
  for(i in 1:length(Twb)){
    # Discrimination of rain, sleet, and snow
    if(Twb[i] > Tmin[i] && Twb[i] < Tmax[i]){
      # sleet
      Pr_sno[i]<-AWS[i,'Rain_TB']*solid_fraction[i];
      Pr_liq[i]<-AWS[i,'Rain_TB']*(1-solid_fraction[i]);
    }
    if(Twb[i] <= Tmin[i]){
      # snow
      Pr_sno[i]<-AWS[i,'Rain_TB']
      Pr_liq[i]<- 0
    }
    if(Twb[i]>= Tmax[i]){
      # rain
      Pr_sno[i]<-0
      Pr_liq[i]<-AWS[i,'Rain_TB']
    }
  }
  
  return(cbind.data.frame(Pr_liq,Pr_sno))
}



###############################################################################
# OPPOSITE Of %IN%
#
'%notin%' <- Negate('%in%')
