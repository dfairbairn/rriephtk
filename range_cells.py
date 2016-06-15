def findRangeCell(lats,lons, myFov=None, elementsPerChunk=5000):
  """ Given input latitudes and longitudes, determine which (if any) 
      SuperDARN range cells contain the latitude/longitude pairs.
  
  **Args**:    
    * **lats** (numpy.ndarray/list): list or array of input latitudes
    * **lons** (numpy.ndarray/list): list or array of input longitudes
    * **[myFov]** (pydarn.radar.radFov.fov): defaults to Saskatoon radar 
                 if none provided.
    * **[elementsPerChunk]** (int): THIS IS EXTREMELY IMPORTANT. SET THIS TO AN 
                          INTEGER THAT DOES NOT USE UP ALL OF YOUR RAM,
                          FOR 16GB of RAM, 50000 USES ~40%.

  **Example**:
    ::
      import range_cells
      beams,gates = range_cells.findRangeCell([60,60],[-100,-90])

  **Note**:
      If myFov is not specified, this code defaults to a Saskatoon FOV
        
  written by A. S. Reimer, 2013-09


  Minutely modified by David Fairbairn, 2016-06

  """
# Given an array of latitude and longitude, this function determines which
# radar cell the point(s) falls within, returning [beam,gate] or [-1,-1] if
# the point is outside the field-of-view.  The routine uses the x and y
# arrays determined by define_beams.  Lat and lon are given in the
# coordinate system set before the last call of this function.
#
# This code has been revised to account for the fact that longitudes can
# flip from -180 to 180 degrees within a radar cell (which confused the old
# code).  So as to cope gracefully with the pole as well, a new approach is
# employed in which everything is converted into a "cartesian/polar" coord
# system, similar to that used for plotting polar diagrams
#
# ASR mod June 14, 2011
# This code has been further revised to handle arrays of lat and lon to be 
# input such that one may process a large number of them efficiently (timewise). 
# since Chapman now has 16GB of RAM we can easily use up some memory (for a 
# short amount of time) to speed up our search. For 483,000 lats and lons, the 
# old way (using for loops) took approximately 6 hours or so, and the new way 
# (vectorization and using the where command) only takes approximately 30 seconds.
#
# Additional modifications (ASR, June 14, 2011)
# ERROR CODES, a -1 means that the latitude and longitude were outside of the 
# field of view and a -2 means that the latitude and longitude arrays input were
# not of the same length.
# 
# This function also now explicitly calls define_beams to ensure that it has the
# correct beam data to work with. The old way was to trust that the user called
# define_beams before calling this function (this violated the principle of 
# least astonishment).
#
# The function has been limited to unsing ~40% of Chapman's RAM (~6.4GB). One may
# limit the memory usage by altering a constant called elementsPerChunk. This breaks 
# the input lat and lon vectors up into elementsPerChunk sized chunks for processing.
#

  from davitpy import pydarn
  import numpy as np
  from datetime import datetime
  import logging

  assert(np.shape(lats) == np.shape(lons)),"Input latitude/longitude must be vectors of the same length."
  #make sure lats and lons are numpy arrays and longitude is the same format as FOV lons (0-360) 
  lons,lats = (np.array(lons)+360.)%360.0, np.array(lats)
  numLatsLons = np.size(lons)
  shapeLatsLons = np.shape(lats)
  lons=lons.reshape((numLatsLons,))
  lats=lats.reshape((numLatsLons,))
  #first let's get the fov of the radar of interest
  if not myFov:
    site = pydarn.radar.site(code='sas')
    myFov = pydarn.radar.radFov.fov(site=site, altitude=300.0, model='IS', coords='geo')

  #ensure FOV lats and lons are in same format as input
  rLons, rLats = (np.array(myFov.lonFull)+360.)%360.0,np.array(myFov.latFull)

  #Now break input up into chucks if necessary and process one chunk at a time
  chunks=int(np.floor((numLatsLons-1)/elementsPerChunk))+1
  rangeCells = np.nan*np.ones((2,numLatsLons))
  numGates=np.size(myFov.gates)
  numBeams=np.size(myFov.beams)
  numFovElements=numGates*numBeams

  sTime=datetime.now()
  #Beginning of loop for each iteration (or chunk of input vector)
  for i in range(0,chunks):
    startInd=i*elementsPerChunk      #Starting index dependent on current iteration
    #but the end index is not
    if (numLatsLons > elementsPerChunk):
      endInd=(i+1)*elementsPerChunk
      numLatsLons=numLatsLons-elementsPerChunk
      numEl=elementsPerChunk      
    else:
      endInd=i*elementsPerChunk+numLatsLons  
      numEl=numLatsLons
    
    logging.info("Processing elements "+str(startInd+1)+" to "+str(endInd)) # Previously a print statement
    
    gates,beams=np.meshgrid(myFov.gates,myFov.beams)
    beams=np.reshape(beams,(1,numFovElements))*np.ones((numEl,1))
    gates=np.reshape(gates,(1,numFovElements))*np.ones((numEl,1))  
  
    #Now we convert fov locations to input position to "cartesian/polar" coords
    #similar to when one wants to plot points on unp.sing the polar plotting routines,
    #except here we don't have to worry about unp.sing the mlt function because we don't 
    #care what time it is.

    posX= (90-np.abs(lats[startInd:endInd]))*np.sin(np.deg2rad(lons[startInd:endInd]))  #now convert them into the  polar/catesian grid
    posY=-(90-np.abs(lats[startInd:endInd]))*np.cos(np.deg2rad(lons[startInd:endInd]))  #coordinate system

    #Next let's create arrays that store the coordinates of the corners of the fov cells
    #and at the same time we convert them into the polar/cartesian coordinate system.
    #BOTTOM LEFT CORNER
    polarX= (90-np.abs(rLats[0:numBeams,0:numGates]))*np.sin(np.deg2rad(rLons[0:numBeams,0:numGates]))
    polarY=-(90-np.abs(rLats[0:numBeams,0:numGates]))*np.cos(np.deg2rad(rLons[0:numBeams,0:numGates]))
    polarX1=np.reshape(polarX,(1,numFovElements))
    polarY1=np.reshape(polarY,(1,numFovElements))

    #TOP RIGHT CORNER
    polarX= (90-np.abs(rLats[1:numBeams+1,1:numGates+1]))*np.sin(np.deg2rad(rLons[1:numBeams+1,1:numGates+1]))
    polarY=-(90-np.abs(rLats[1:numBeams+1,1:numGates+1]))*np.cos(np.deg2rad(rLons[1:numBeams+1,1:numGates+1]))
    polarX3=np.reshape(polarX,(1,numFovElements))
    polarY3=np.reshape(polarY,(1,numFovElements))

    #TOP LEFT CORNER
    polarX= (90-np.abs(rLats[0:numBeams,1:numGates+1]))*np.sin(np.deg2rad(rLons[0:numBeams,1:numGates+1]))
    polarY=-(90-np.abs(rLats[0:numBeams,1:numGates+1]))*np.cos(np.deg2rad(rLons[0:numBeams,1:numGates+1]))
    polarX4=np.reshape(polarX,(1,numFovElements))
    polarY4=np.reshape(polarY,(1,numFovElements))

    #BOTTOM RIGHT CORNER
    polarX= (90-np.abs(rLats[1:numBeams+1,0:numGates]))*np.sin(np.deg2rad(rLons[1:numBeams+1,0:numGates]))
    polarY=-(90-np.abs(rLats[1:numBeams+1,0:numGates]))*np.cos(np.deg2rad(rLons[1:numBeams+1,0:numGates]))
    polarX2=np.reshape(polarX,(1,numFovElements))
    polarY2=np.reshape(polarY,(1,numFovElements))

    #Now, due to how the where command works, we need to convert the arrays that we just
    #made into matricies. This ensures that we compare every input lat and lon with every 
    #corner coordinate of the fov of the radar. See the where command documentation for 
    #more details. 

    logging.info("Setting up comparison matricies...") # Previously a print statement

    catX2=polarX2*np.ones((numEl,1))
    catY2=polarY2*np.ones((numEl,1))
    catX1=polarX1*np.ones((numEl,1))
    catY1=polarY1*np.ones((numEl,1))
    catX3=polarX3*np.ones((numEl,1))
    catY3=polarY3*np.ones((numEl,1))
    catX4=polarX4*np.ones((numEl,1))
    catY4=polarY4*np.ones((numEl,1))
    catX0=np.ones((1,numFovElements))*np.reshape(posX,(numEl,1))
    catY0=np.ones((1,numFovElements))*np.reshape(posY,(numEl,1))

    #To find the overlaps, now we generate the line equation that is used to draw the cell
    #this requires us to assume that the cells are made of straight lines (not valid on a 
    #sphere, but we are working in a cartesian/polar coordinate system so it is valid).
    #Then we test to make sure that the latitude and longitude points are bounded by these
    #line equations.

    logging.info("Finding overlaps...") # Previously a print statement

    m1=(catY2-catY1)/(catX2-catX1)    #first we need to find the slopes of the lines
    m2=(catY3-catY2)/(catX3-catX2)    #by using the corners of the cell
    m3=(catY4-catY3)/(catX4-catX3)
    m4=(catY1-catY4)/(catX1-catX4)
      
    b1=catY2-m1*catX2        #next we calculate the y-intercept of the lines
    b2=catY3-m2*catX3
    b3=catY4-m3*catX4
    b4=catY1-m4*catX1
      
    x1=(catY0-b1)/m1        #now let's calculate the x value that would result
    x2=(catY0-b2)/m2        #if we plug in the input y value (input y value was 
    x3=(catY0-b3)/m3        #calculated from the input lat and lon)
    x4=(catY0-b4)/m4
      
      
    #SPECIAL CASE (slope was infinite, so x1,x2,x3,x4 ended up being infinte too)
    ind=np.where(x1 == np.nan)[0]    #So then we must make our x value generated from the input y value
    if (np.size(ind) > 0):              #equal to the x value of one of the corners
      x1[ind]=catX2[ind]
    ind=np.where(x2 == np.nan)[0]
    if (np.size(ind) > 0):
      x2[ind]=catX3[ind]
    ind=np.where(x3 == np.nan)[0]
    if (np.size(ind) > 0):
      x3[ind]=catX4[ind]
    ind=np.where(x4 == np.nan)[0]
    if (np.size(ind) > 0):
      x4[ind]=catX1[ind]
  
    #Now, the magic of the where command comes in.
    #We must now check to see if the input x and y values are bounded by the lines of the cell
      
    ind=np.where((((np.abs(catX1-catX2) > np.abs(catX1-x1)) & (np.abs(catX1-catX2) > np.abs(x1-catX2))) & \
                  ((np.abs(catX2-catX3) > np.abs(catX2-x2)) & (np.abs(catX2-catX3) > np.abs(x2-catX3))) & \
                  ((np.abs(x1-x2) > np.abs(x1-catX0)) & (np.abs(x1-x2) > np.abs(x2-catX0)))) \
                  | \
                 (((np.abs(catX1-catX2) > np.abs(catX1-x1)) & (np.abs(catX1-catX2) > np.abs(x1-catX2))) & \
                  ((np.abs(catX3-catX4) > np.abs(catX3-x3)) & (np.abs(catX3-catX4) > np.abs(x3-catX4))) & \
                  ((np.abs(x1-x3) > np.abs(x1-catX0)) & (np.abs(x1-x3) > np.abs(x3-catX0)))) \
                  | \
                 (((np.abs(catX1-catX2) > np.abs(catX1-x1)) & (np.abs(catX1-catX2) > np.abs(x1-catX2))) & \
                  ((np.abs(catX4-catX1) > np.abs(catX4-x4)) & (np.abs(catX4-catX1) > np.abs(x4-catX1))) & \
                  ((np.abs(x1-x4) > np.abs(x1-catX0)) & (np.abs(x1-x4) > np.abs(x4-catX0)))) \
                  | \
                 (((np.abs(catX2-catX3) > np.abs(catX2-x2)) & (np.abs(catX2-catX3) > np.abs(x2-catX3))) & \
                  ((np.abs(catX3-catX4) > np.abs(catX3-x3)) & (np.abs(catX3-catX4) > np.abs(x3-catX4))) & \
                  ((np.abs(x2-x3) > np.abs(x2-catX0)) & (np.abs(x2-x3) > np.abs(x3-catX0)))) \
                  | \
                 (((np.abs(catX2-catX3) > np.abs(catX2-x2)) & (np.abs(catX2-catX3) > np.abs(x2-catX3))) & \
                  ((np.abs(catX4-catX1) > np.abs(catX4-x4)) & (np.abs(catX4-catX1) > np.abs(x4-catX1))) & \
                  ((np.abs(x2-x4) > np.abs(x2-catX0)) & (np.abs(x2-x4) > np.abs(x4-catX0)))) \
                  | \
                 (((np.abs(catX3-catX4) > np.abs(catX3-x3)) & (np.abs(catX3-catX4) > np.abs(x3-catX4))) & \
                  ((np.abs(catX4-catX1) > np.abs(catX4-x4)) & (np.abs(catX4-catX1) > np.abs(x4-catX1))) & \
                  ((np.abs(x3-x4) > np.abs(x3-catX0)) & (np.abs(x3-x4) > np.abs(x4-catX0)))))

    if (np.size(ind[0]) > 0): #if np.size(ind) > 0 then we found some overlaps! So now put them into the output vector.
      rangeCells[0, startInd+ind[0]]=beams[ind[0],ind[1]]
      rangeCells[1, startInd+ind[0]]=gates[ind[0],ind[1]]
  fTime=datetime.now()-sTime

  logging.info("Finished in "+ str(fTime)) # Previously a print statement


# So we had to do a bit more thinking, and a bit more programming 
# (the code is longer than only using for loops) but the result
# is that the function is about a thousand times faster!

  return np.reshape(rangeCells[0,:],shapeLatsLons), np.reshape(rangeCells[1,:],shapeLatsLons),
