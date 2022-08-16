### SLATED FOR DELETION; PROBABLY NO LONGER NEEDED, BUT FIRST SOME FUNCTIONS SHOULD BE MOVED TO OTHER SCRIPTS
# ----------------------------------------------------------------------------------------
# libScuFx.py
# Version:  ArcGIS 10.3.1 / Python 2.7.8
# Creation Date: 2017-08-29
# Last Edit: 2018-11-05
# Creator(s):  Kirsten R. Hazler

# Summary:
# A library of functions for delineating and prioritizing Stream Conservation Units (SCUs) for conservation.

# Usage Tips:
# 

# Dependencies:
# 

# Syntax:  
# 
# ----------------------------------------------------------------------------------------

# Import modules
import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
import libConSiteFx
from libConSiteFx import *
import os, sys, datetime, traceback, gc

###DELETE THIS FUNCTION AFTER SUCCESSFULLY MOVED TO CreateSCU.py
def delinFlowDistBuff(in_Feats, fld_ID, in_FlowDir, out_Feats, maxDist, dilDist = 0, out_Scratch = 'in_memory'):
   '''Delineates buffers based on flow distance down to features (rather than straight distance)'''
   # Get cell size and output spatial reference from in_FlowDir
   cellSize = (arcpy.GetRasterProperties_management(in_FlowDir, "CELLSIZEX")).getOutput(0)
   srRast = arcpy.Describe(in_FlowDir).spatialReference
   linUnit = srRast.linearUnitName
   printMsg('Cell size of flow direction raster is %s %ss' %(cellSize, linUnit))
   printMsg('Flow modeling is strongly dependent on cell size.')

   # Set environment setting and other variables
   arcpy.env.overwriteOutput = True
   arcpy.env.snapRaster = in_FlowDir
   procDist = 2*maxDist
   # dist = float(cellSize)
   # dilDist = "%s %ss" % (str(dist), linUnit) # decided to make this a user input

   # Check if input features and input flow direction have same spatial reference.
   # If so, just make a copy. If not, reproject features to match raster.
   srFeats = arcpy.Describe(in_Feats).spatialReference
   if srFeats.Name == srRast.Name:
      printMsg('Coordinate systems for features and raster are the same. Copying...')
      arcpy.CopyFeatures_management (in_Feats, out_Feats)
   else:
      printMsg('Reprojecting features to match raster...')
      # Check if geographic transformation is needed, and handle accordingly.
      if srFeats.GCS.Name == srRast.GCS.Name:
         geoTrans = ""
         printMsg('No geographic transformation needed...')
      else:
         transList = arcpy.ListTransformations(srFeats,srRast)
         geoTrans = transList[0]
      arcpy.Project_management (in_Feats, out_Feats, srRast, geoTrans)

   # Count features and report
   numFeats = countFeatures(out_Feats)
   printMsg('There are %s features to process.' % numFeats)
      
   # Create an empty list to store IDs of features that fail to get processed
   myFailList = []

   # Set up processing cursor and loop
   flags = [] # Initialize empty list to keep track of suspects
   cursor = arcpy.da.UpdateCursor(out_Feats, [fld_ID, "SHAPE@"])
   counter = 1
   for row in cursor:
      trashList = [] # Empty list for trash collection
      try:
         # Extract the unique ID and geometry object
         myID = row[0]
         myShape = row[1]

         printMsg('Working on feature %s with ID %s' % (counter, str(myID)))

         # Process:  Select (Analysis)
         # Create a temporary feature class including only the current feature
         selQry = "%s = %s" % (fld_ID, str(myID))
         tmpFeat = out_Scratch + os.sep + 'tmpFeat'
         trashList.append(tmpFeat)
         arcpy.Select_analysis (out_Feats, tmpFeat, selQry)

         # Convert feature to raster
         printMsg('Converting feature to raster...')
         srcRast = out_Scratch + os.sep + 'srcRast'
         trashList.append(srcRast)
         arcpy.PolygonToRaster_conversion (tmpFeat, fld_ID, srcRast, "MAXIMUM_COMBINED_AREA", fld_ID, cellSize)
         
         # Clip flow direction raster to processing buffer
         procBuff = out_Scratch + os.sep + 'procBuff'
         trashList.append(procBuff)
         printMsg('Buffering feature to set maximum processing distance')
         arcpy.Buffer_analysis (tmpFeat, procBuff, procDist, "", "", "ALL", "")
         myExtent = str(arcpy.Describe(procBuff).extent).replace(" NaN", "")
         printMsg('Extent: %s' %myExtent)
         clp_FlowDir = out_Scratch + os.sep + 'clp_FlowDir'
         trashList.append(clp_FlowDir)
         printMsg('Clipping flow direction raster to processing buffer')
         arcpy.Clip_management (in_FlowDir, myExtent, clp_FlowDir, procBuff, "", "ClippingGeometry")
         arcpy.env.extent = clp_FlowDir

         # Burn SCU feature into flow direction raster as sink
         printMsg('Creating sink from feature...')
         snk_FlowDir = Con(IsNull(srcRast),clp_FlowDir)
         trashList.append(snk_FlowDir)
         snk_FlowDir.save(out_Scratch + os.sep + 'snk_FlowDir')
         
         # Calculate flow distance down to sink
         printMsg('Calculating flow distance to feature...')
         FlowDist = FlowLength (snk_FlowDir, "DOWNSTREAM")
         FlowDist.save(out_Scratch + os.sep + 'FlowDist')
         trashList.append(FlowDist)
         
         # Clip flow distance raster to the maximum distance buffer
         clipBuff = out_Scratch + os.sep + 'clipBuff'
         trashList.append(clipBuff)
         arcpy.Buffer_analysis (tmpFeat, clipBuff, maxDist, "", "", "ALL", "")
         myExtent = str(arcpy.Describe(clipBuff).extent).replace(" NaN", "")
         #printMsg('Extent: %s' %myExtent)
         clp_FlowDist = out_Scratch + os.sep + 'clp_FlowDist'
         trashList.append(clp_FlowDist)
         printMsg('Clipping flow distance raster to maximum distance buffer')
         arcpy.Clip_management (FlowDist, myExtent, clp_FlowDist, clipBuff, "", "ClippingGeometry")
         arcpy.env.extent = clp_FlowDist
         
         # Make a binary raster based on flow distance
         printMsg('Creating binary raster from flow distance...')
         binRast = Con((IsNull(clp_FlowDist) == 1),
                  (Con((IsNull(srcRast)== 0),1,0)),
                  (Con((Raster(clp_FlowDist) <= maxDist),1,0)))
         binRast.save(out_Scratch + os.sep + 'binRast')
         trashList.append(binRast)
         printMsg('Boundary cleaning...')
         cleanRast = BoundaryClean (binRast, 'NO_SORT', 'TWO_WAY')
         cleanRast.save(out_Scratch + os.sep + 'cleanRast')
         trashList.append(cleanRast)
         printMsg('Setting zeros to nulls...')
         prePoly = SetNull (cleanRast, 1, 'Value = 0')
         prePoly.save(out_Scratch + os.sep + 'prePoly')
         trashList.append(prePoly)

         # Convert raster to polygon
         printMsg('Converting flow distance raster to polygon...')
         finPoly = out_Scratch + os.sep + 'finPoly'
         trashList.append(finPoly)
         arcpy.RasterToPolygon_conversion (prePoly, finPoly, "NO_SIMPLIFY")

         # If user specifies, coalesce to smooth
         if dilDist == 0:
            printMsg('Final shape will not be smoothed.')
            coalPoly = finPoly
         else:
            printMsg('Smoothing final shape...')
            coalPoly = out_Scratch + os.sep + 'coalPoly'
            Coalesce(finPoly, dilDist, coalPoly, 'in_memory')
         trashList.append(coalPoly)
         
         # Check the number of features at this point. 
         # It should be just one. If more, the output is likely bad and should be flagged.
         count = countFeatures(coalPoly)
         if count > 1:
            printWrng('Output is suspect for feature %s' % str(myID))
            flags.append(myID)
            # Dissolve to create a multipart feature so at least we can finish the job.
            multiPoly = out_Scratch + os.sep + 'multiPoly'
            arcpy.Dissolve_management (coalPoly, multiPoly, "", "", "MULTI_PART")
            coalPoly = multiPoly
         
         # Use the flow distance buffer geometry as the final shape
         myFinalShape = arcpy.SearchCursor(coalPoly).next().Shape

         # Update the feature with its final shape
         row[1] = myFinalShape
         cursor.updateRow(row)
         del row 

         printMsg('Finished processing feature %s' %str(myID))
         
         # Reset extent, because Arc is stupid.
         arcpy.env.extent = "MAXOF"
         
         # Throw out trash after every cycle and compact the scratch GDB periodically. 
         # Grasping at straws here to avoid failure processing large datasets.
         garbagePickup(trashList)        
         for item in trashList:
            del item
         if counter%25 == 0:
            printMsg('Compacting scratch geodatabase...')
            arcpy.Compact_management (out_Scratch)
         
         # Update counter
         counter += 1

      except:
         # Add failure message and append failed feature ID to list
         printMsg("\nFailed to fully process feature " + str(myID))
         myFailList.append(myID)

         # Error handling code swiped from "A Python Primer for ArcGIS"
         tb = sys.exc_info()[2]
         tbinfo = traceback.format_tb(tb)[0]
         pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n " + str(sys.exc_info()[1])
         msgs = "ARCPY ERRORS:\n" + arcpy.GetMessages(2) + "\n"

         printWrng(msgs)
         printWrng(pymsg)
         printMsg(arcpy.GetMessages(1))

         # Add status message
         printMsg("\nMoving on to the next feature.  Note that the output will be incomplete.")
   
   if len(flags) > 0:
      printWrng('These features may be incorrect: %s' % str(flags))
   if len(myFailList) > 0:
      printWrng('These features failed to process: %s' % str(myFailList))
   return out_Feats

def getZonalStats(in_Polys, in_Raster, fld_ID, fld_Stats, type_Stats, out_Polys, out_Scratch = 'in_memory'):
   '''Attaches zonal statistics from a raster to polygons'''
   # Environment settings
   arcpy.env.overwriteOutput = True
   cellSize = (arcpy.GetRasterProperties_management(in_Raster, "CELLSIZEX")).getOutput(0)   
   arcpy.env.cellSize = cellSize
   printMsg('Cell Size is %s map units' % str(cellSize))
   arcpy.env.snapRaster = in_Raster
   
   # Check if input polygons and raster have same spatial reference.
   # If so, all good. If not, reproject features to match raster.
   srPolys = arcpy.Describe(in_Polys).spatialReference
   srRast = arcpy.Describe(in_Raster).spatialReference
   if srPolys.Name == srRast.Name:
      printMsg('Coordinate systems for features and raster are the same. Copying...')
      arcpy.CopyFeatures_management (in_Polys, out_Polys)
   else:
      printMsg('Reprojecting features to match raster...')
      # Check if geographic transformation is needed, and handle accordingly.
      if srPolys.GCS.Name == srRast.GCS.Name:
         geoTrans = ""
         printMsg('No geographic transformation needed...')
      else:
         transList = arcpy.ListTransformations(srPolys,srRast)
         geoTrans = transList[0]
      arcpy.Project_management (in_Polys, out_Polys, srRast, geoTrans)
      
   # Set up some variables used in loop below
   tmpFeat = out_Scratch + os.sep + 'tmpFeat' # temp feature class 
   tmpTab = out_Scratch + os.sep + 'tmpTab' # temp table 
   mstrTab = out_Scratch + os.sep + 'mstrTab' # master table
   for t in [tmpTab, mstrTab]:
      if arcpy.Exists(t):
         arcpy.Delete_management(t)
   myFailList = [] # empty list to keep track of failures
   
   # Count features and report
   numFeats = countFeatures(out_Polys)
   printMsg('There are %s features to process.' % numFeats)
   
   # Set up processing cursor and loop through polygons to get zonal stats. 
   # This is done in a loop instead of all at once to avoid problems if polys overlap.
   cursor = arcpy.da.SearchCursor(out_Polys, '%s' % fld_ID)
   counter = 1
   for row in cursor:
      try:
         myID = row[0]
         printMsg('Working on zonal stats for feature %s with ID %s' % (counter, str(myID)))
         
         # Make temporary feature class
         selQry = "%s = %s" % (fld_ID, str(myID))
         arcpy.Select_analysis (out_Polys, tmpFeat, selQry)
         
         # Run zonal stats and append to master table
         outTab = ZonalStatisticsAsTable (tmpFeat, fld_ID, in_Raster, tmpTab, 'DATA', 'ALL')
         if arcpy.Exists(mstrTab):
            arcpy.Append_management(tmpTab, mstrTab, 'NO_TEST')
         else:
            arcpy.Copy_management(tmpTab, mstrTab)
         
         # Updates
         del row, outTab
         counter += 1
                  
         # Reset extent, because Arc is stupid.
         arcpy.env.extent = "MAXOF"
         
      except:
      # Add failure message and append failed feature ID to list
         printMsg("\nFailed to fully process feature " + str(myID))
         myFailList.append(myID)

         # Error handling code swiped from "A Python Primer for ArcGIS"
         tb = sys.exc_info()[2]
         tbinfo = traceback.format_tb(tb)[0]
         pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n " + str(sys.exc_info()[1])
         msgs = "ARCPY ERRORS:\n" + arcpy.GetMessages(2) + "\n"

         printWrng(msgs)
         printWrng(pymsg)
         printMsg(arcpy.GetMessages(1))

         # Add status message
         printMsg("\nMoving on to the next feature.  Note that the output will be incomplete.")
   del cursor
   
   # For each stats type, manipulate the fields
   printMsg('Appending stats fields to output features')
   polyFldNames = [field.name for field in arcpy.ListFields(out_Polys)]
   printMsg(polyFldNames)
   for t in type_Stats:
      fldName = fld_Stats + '_' + t
      if fldName in polyFldNames:
         arcpy.DeleteField_management(out_Polys, fldName)
      arcpy.AddField_management (mstrTab, fldName, 'FLOAT')
      expression = '!%s!' % t
      try:
         arcpy.CalculateField_management(mstrTab, fldName, expression, 'PYTHON')
         arcpy.JoinField_management(out_Polys, fld_ID, mstrTab, fld_ID, fldName)
         # I put this in a try block b/c certain stats, such as MEDIAN, are apparently only available for integer rasters.
      except:
         printWrng('Unable to add %s' %t)

         # Error handling code swiped from "A Python Primer for ArcGIS"
         tb = sys.exc_info()[2]
         tbinfo = traceback.format_tb(tb)[0]
         pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n " + str(sys.exc_info()[1])
         msgs = "ARCPY ERRORS:\n" + arcpy.GetMessages(2) + "\n"

         printWrng(msgs)
         printWrng(pymsg)
         printMsg(arcpy.GetMessages(1))

   return out_Polys
   
def getLandscapeScore(in_Feats, fld_ForWet, fld_ImpSur):
   '''Scores features based on forest/wetland and impervious surface cover'''
   
   # If needed, add fields to hold the scores
   printMsg('Adding fields')
   fldNames = [field.name for field in arcpy.ListFields(in_Feats)]
   for f in ['ForWet_score', 'ImpSur_score', 'Landscape_score']:
      if f not in fldNames:
         arcpy.AddField_management (in_Feats, f, 'FLOAT')
      
   # Get forest/wetland score
   printMsg('Calculating forest/wetland score')
   codeblock = '''def forwetScore(forwet):
      fMax = 0.85
      fMin = 0.15
      if forwet <= fMin:
         score = 0
      elif forwet >= fMax:
         score = 100
      else:
         score = 100*(forwet - fMin)/(fMax - fMin)
      return score'''
   expression = 'forwetScore(!%s!)' % fld_ForWet
   arcpy.CalculateField_management(in_Feats, 'ForWet_score', expression, "PYTHON", codeblock)
         
   # Get impervious surface score
   printMsg('Calculating impervious surface score')
   codeblock = '''def impsurScore(impsur):
      iMax = 0.25
      iMin = 0.05
      if impsur <= iMin:
         score = 100
      elif impsur >= iMax:
         score = 0
      else:
         score = 100*(iMax - impsur)/(iMax - iMin)
      return score''' 
   expression = 'impsurScore(!%s!)' % fld_ImpSur
   arcpy.CalculateField_management(in_Feats, 'ImpSur_score', expression, "PYTHON", codeblock)
     
   # Get landscape score
   printMsg('Calculating landscape score')
   expression = r'(!ForWet_score! + !ImpSur_score!)/2' 
   arcpy.CalculateField_management(in_Feats, 'Landscape_score', expression, "PYTHON")
   
   return in_Feats

def prioritizeSCUs(in_Feats, fld_ID, fld_BRANK, lo_BRANK, in_Integrity, lo_Integrity, in_ConsPriority, in_Vulnerability, out_Feats, out_Scratch = 'in_memory'):
   '''Prioritizes Stream Conservation Units (SCUs) for conservation, based on biodiversity rank (BRANK), watershed integrity and conservation priority (from ConservationVision Watershed Model), and vulnerability (from ConservationVision Development Vulnerability Model)'''
   
   # Environment settings
   arcpy.env.overwriteOutput = True
   cellSize = (arcpy.GetRasterProperties_management(in_Integrity, "CELLSIZEX")).getOutput(0)   
   arcpy.env.cellSize = cellSize
   printMsg('Cell Size is %s map units' % str(cellSize))
   arcpy.env.snapRaster = in_Integrity
   
   # Step 1: First cut: Create subset of buffered SCUs ranked lo_BRANK or better
   selQry = "%s <= '%s'" % (fld_BRANK, lo_BRANK)
   arcpy.Select_analysis (in_Feats, out_Feats, selQry)
   
   # Step 2: For each buffered SCU in subset, get zonal stats for Watershed Integrity, Conservation Priority and Vulnerability. Do this in loop in case of buffer overlap.
   
   # Add fields to hold the stats
   for f in ['INTEG', 'CPRIOR_MEAN', 'CPRIOR_MAX', 'VULN', 'SCORE']:
      arcpy.AddField_management (out_Feats, f, 'FLOAT')
      
   # Set up some variables used in loop
   tmpFeat = out_Scratch + os.sep + 'tmpFeat' # temp feature class 
   tmpTab = out_Scratch + os.sep + 'tmpTab' # temp table 
   if arcpy.Exists(tmpTab):
      arcpy.Delete_management(tmpTab)
   myFailList = [] # empty list to keep track of failures
   
   # Count features and report
   numFeats = countFeatures(out_Feats)
   printMsg('There are %s features to process.' % numFeats)
   
   # Set up processing cursor and loop
   fields = [fld_ID, 'INTEG', 'CPRIOR_MEAN', 'CPRIOR_MAX', 'VULN']
   cursor = arcpy.da.UpdateCursor(out_Feats, fields)
   counter = 1
   for row in cursor:
      try:
         myID = row[0]
         
         printMsg('Working on zonal stats for feature %s with ID %s' % (counter, str(myID)))
         
         # Make temporary feature class
         selQry = "%s = %s" % (fld_ID, str(myID))
         arcpy.Select_analysis (out_Feats, tmpFeat, selQry)
         
         # Run zonal stats, extract values, and update fields
         outTab = ZonalStatisticsAsTable (tmpFeat, fld_ID, in_Integrity, tmpTab, 'DATA', 'MEAN')
         val = arcpy.SearchCursor(tmpTab).next().MEAN
         printMsg('The watershed integrity mean is %s' % str(val))
         row[1] = val
         
         outTab = ZonalStatisticsAsTable (tmpFeat, fld_ID, in_ConsPriority, tmpTab, 'DATA', 'ALL')
         val = arcpy.SearchCursor(tmpTab).next().MEAN
         printMsg('The conservation priority mean is %s' % str(val))
         row[2] = val
         val = arcpy.SearchCursor(tmpTab).next().MAX
         printMsg('The conservation priority maximum is %s' % str(val))
         row[3] = val
         
         outTab = ZonalStatisticsAsTable (tmpFeat, fld_ID, in_Vulnerability, tmpTab, 'DATA', 'MEAN')
         val = arcpy.SearchCursor(tmpTab).next().MEAN
         printMsg('The development vulnerability mean is %s' % str(val))
         row[4] = val
         
         # Updates
         cursor.updateRow(row)
         del row, val, outTab
         counter += 1
                  
         # Reset extent, because Arc is stupid.
         arcpy.env.extent = "MAXOF"
         
      except:
      # Add failure message and append failed feature ID to list
         printMsg("\nFailed to fully process feature " + str(myID))
         myFailList.append(myID)

         # Error handling code swiped from "A Python Primer for ArcGIS"
         tb = sys.exc_info()[2]
         tbinfo = traceback.format_tb(tb)[0]
         pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n " + str(sys.exc_info()[1])
         msgs = "ARCPY ERRORS:\n" + arcpy.GetMessages(2) + "\n"

         printWrng(msgs)
         printWrng(pymsg)
         printMsg(arcpy.GetMessages(1))

         # Add status message
         printMsg("\nMoving on to the next feature.  Note that the output will be incomplete.")

   # Step 3: Score catchments based on BRANK, Watershed Integrity, Conservation Priority, and Vulnerability
   codeblock = '''def ConScore(brank, integ, cprior_mean, cprior_max, vuln):
      bDict = {'B1': 5, 'B2': 0, 'B3': -5, 'B4': -10, 'B5': -15}
      b = bDict[brank]
      if integ >= %s:
         i = 1
      else:
         i = 0
      if vuln >= 80:
         v = 0
      elif vuln <= 20:
         v = 20
      else:
         v = -vuln/3 + 26.67
      cprior = (cprior_mean + cprior_max)/2
      rawScore = i*((integ + cprior)/2 + b - v)
      if rawScore > 100:
         score = 100
      elif rawScore < 0:
         score = 0
      else: 
         score = rawScore
      return score''' % lo_Integrity
   expression = 'ConScore(!%s!, !INTEG!, !CPRIOR_MEAN!, !CPRIOR_MAX!, !VULN!)' % fld_BRANK
   arcpy.CalculateField_management(out_Feats, 'SCORE', expression, "PYTHON", codeblock)
   
   if len(myFailList) > 0:
      printWrng('These features failed to process: %s' % str(myFailList))
      
   return out_Feats
   
# Use the main function below to run the delinFlowDistBuff function directly from Python IDE or command line with hard-coded variables
def main():
   in_Feats = r'C:\Users\xch43889\Documents\ArcGIS\Default.gdb\scuBaseline_Dissolve'
   fld_ID = 'ID'
   in_FlowDir = r'H:\Backups\DCR_Work_DellD\GIS_Data_VA_proc\Finalized\NHDPlus_Virginia.gdb\fdir_VA'
   out_Feats = r'C:\Users\xch43889\Documents\ArcGIS\Default.gdb\scuBaseline_FlowBuff1k'
   maxDist = 1000
   dilDist = 0
   out_Scratch = r'C:\Users\xch43889\Documents\ArcGIS\scratch.gdb'
   # End of user input

   delinFlowDistBuff(in_Feats, fld_ID, in_FlowDir, out_Feats, maxDist, dilDist, out_Scratch)
   
# # Use the main function below to run the getZonalStats function
# def main():
   # # Set up hard-coded variables
   # in_Polys = r'C:\Users\xch43889\Documents\Working\SCU_prioritization\SCU_work_20170903.gdb\scuFlowBuff250_examples'
   # in_Raster = r'C:\Users\xch43889\Documents\Working\SCU_prioritization\ConsPrior.tif'
   # fld_ID = 'lngID'
   # fld_Stats = 'ConsPrior'
   # type_Stats = ['MIN', 'MAX', 'MEAN', 'MEDIAN']
   # out_Polys = r'C:\Users\xch43889\Documents\Working\SCU_prioritization\scratch.gdb\test'
   # out_Scratch = r'C:\Users\xch43889\Documents\Working\SCU_prioritization\scratch.gdb'
   # # End of user input
   
   # getZonalStats(in_Polys, in_Raster, fld_ID, fld_Stats, type_Stats, out_Polys, out_Scratch)

if __name__ == '__main__':
   main()
