
####################################################################
def ProjectImagingMesh(work_dir,data_set):
  import vtk
  import vtk.util.numpy_support as vtkNumPy 
  import numpy
  import scipy.io as scipyio
  # echo vtk version info
  print "using vtk version", vtk.vtkVersion.GetVTKVersion()
  # FIXME  notice that order of operations is IMPORTANT
  # FIXME   translation followed by rotation will give different results
  # FIXME   than rotation followed by translation
  # FIXME  Translate -> RotateZ -> RotateY -> RotateX -> Scale seems to be the order of paraview
  vtkReader = vtk.vtkDataSetReader() 
  # offset averaging +/- in both directions about center
  listoffset = [-2,-1,0,1,2]
  if   (data_set == "dog1"):
    # should be in meters
    #AffineTransform.Translate(0.206,0.289,-0.0215)
    #AffineTransform.RotateZ( 0.0 )
    #AffineTransform.RotateY( 0.0 )
    #AffineTransform.RotateX( 90.0 )
    Translate = (0.293,0.0634,0.12345)
    RotateZ   =  -90.0 
    RotateY   =  0.0 
    RotateX   =  90.0 
    Scale     = [1.,1.,1.]
    ImageFileTemplate  = '/FUS4/data2/BioTex/Brain_PhaseII/060626_BrainDog1_Treat/ProcessedTMAPData/BioTexCanines/dog1dtmap.0000.vtk'
  elif (data_set == "dog2"):
    # should be in meters
    Translate = (0.2335,0.235,-0.0335)
    RotateZ   =  6.0  
    RotateY   =  0.0  
    RotateX   =  90.0 
    Scale     = [1.,1.,1.]
    ImageFileTemplate  = '/data/sjfahrenholtz/mdacc/Dog/dog2_98000/mrivis/temperature.0000.vtk' 
    ImageFileTemplate  = '/FUS4/data2/BioTex/Brain_PhaseII/060626_BrainDog1_Treat/ProcessedTMAPData/BioTexCanines/dog2dtmap.0000.vtk'
  elif (data_set == "dog3"):
    # should be in meters
    Translate  = (0.209,0.2625,-0.0247)
    RotateZ    =  -8.0 
    RotateY    =  0.0 
    RotateX    =  90.0 
    Scale      =  [1.,1.,1.]
    ImageFileTemplate  = '/data/sjfahrenholtz/mdacc/Dog/dog3_98000/mrivis/temperature.0000.vtk' 
    ImageFileTemplate  = '/FUS4/data2/BioTex/Brain_PhaseII/060626_BrainDog1_Treat/ProcessedTMAPData/BioTexCanines/dog3dtmap.0000.vtk'
  elif (data_set == "dog4"):
    # should be in meters
    Translate = (0.195,0.245,-0.0715)
    RotateZ   =  -15.0 
    RotateY   =  0.0 
    RotateX   =  90.0 
    Scale     = [1.,1.,1.]
    ImageFileTemplate  = '/FUS4/data2/BioTex/Brain_PhaseII/060626_BrainDog1_Treat/ProcessedTMAPData/BioTexCanines/dog4dtmap.0000.vtk'
  elif (data_set == "human0"):
    # should be in meters
    vtkReader = vtk.vtkXMLImageDataReader() 
    Translate = (0.051,0.080,0.0509)
    RotateZ   =  29.0 
    RotateY   =  86.0 
    RotateX   =  0.0 
    Scale     = [1.,1.,1.]
    ImageFileTemplate  = '/data/fuentes/biotex/090318_751642_treat/Processed/imaging/temperature.0000.vti' 
  elif (data_set == "agar0"):
    # should be in meters
    Translate = (0.0,0.13,-0.095)
    RotateZ   =  0.0 
    RotateY   =  180.0 
    RotateX   =  0.0 
    Scale     = [1.,1.,1.]
    ImageFileTemplate  = '/data/jyung/MotionCorrection/motion_phantom_sept11/tmap.0000.vtk' 
  else:
    raise RuntimeError("\n\n unknown case... ")
  
  # read imaging data geometry that will be used to project FEM data onto
  vtkReader.SetFileName( ImageFileTemplate ) 
  vtkReader.Update()
  templateImage = vtkReader.GetOutput()
  dimensions = templateImage.GetDimensions()
  spacing = templateImage.GetSpacing()
  origin  = templateImage.GetOrigin()
  print spacing, origin, dimensions
  #fem.SetImagingDimensions( dimensions ,origin,spacing) 

  #setup to interpolate at 5 points across axial dimension 
  TransformList = []
  naverage = len(listoffset)
  subdistance = spacing[2] / (naverage-1)
  print "subspacing distance = ", subdistance
  for idtransform in listoffset: 
    AffineTransform = vtk.vtkTransform()
    AffineTransform.Translate( Translate[0],Translate[1],
                               Translate[2] + subdistance*idtransform )
    AffineTransform.RotateZ( RotateZ )
    AffineTransform.RotateY( RotateY )
    AffineTransform.RotateX( RotateX )
    AffineTransform.Scale(   Scale   )
    TransformList.append( AffineTransform )
    #laserTip         =  AffineTransform.TransformPoint(  laserTip )
    #laserOrientation =  AffineTransform.TransformVector( laserOrientation )

  # Interpolate FEM onto imaging data structures
  vtkExodusIIReader = vtk.vtkExodusIIReader()
  #vtkExodusIIReader.SetFileName( "%s/fem_stats.e" % work_dir )
  #vtkExodusIIReader.SetPointResultArrayStatus("Mean0",1)
  #vtkExodusIIReader.SetPointResultArrayStatus("StdDev0",1)
  #Timesteps  = 120
  for ii in range(0,120):
	  vtkExodusIIReader.SetFileName( "%s/fem_stats.%04d.e" % (work_dir,ii) )
	  #vtkExodusIIReader.SetPointResultArrayStatus("Vard0Mean",1)
	  #vtkExodusIIReader.SetPointResultArrayStatus("Vard0Kurt",1)
	  vtkExodusIIReader.Update()	  
	  numberofresultarrays = vtkExodusIIReader.GetNumberOfPointResultArrays()
	  print numberofresultarrays
	  for resultarrayindex in range(numberofresultarrays):
		resultarrayname = vtkExodusIIReader.GetPointResultArrayName(resultarrayindex)
		vtkExodusIIReader.SetPointResultArrayStatus( "%s" % (resultarrayname),1)
		print resultarrayname
	  vtkExodusIIReader.Update()
	  ntime  = vtkExodusIIReader.GetNumberOfTimeSteps()
	  #for timeID in range(69,70):
	  for timeID in range(ntime):
	    vtkExodusIIReader.SetTimeStep(timeID) 
	    vtkExodusIIReader.Update()
	    
	    # reflect
	    vtkReflectX = vtk.vtkReflectionFilter()
	    vtkReflectX.SetPlaneToXMin()
	    vtkReflectX.SetInput( vtkExodusIIReader.GetOutput() )
	    vtkReflectX.Update()

	    # reflect
	    vtkReflectY = vtk.vtkReflectionFilter()
	    vtkReflectY.SetPlaneToYMax()
	    vtkReflectY.SetInput( vtkReflectX.GetOutput() )
	    vtkReflectY.Update()

	    # apply the average of the transform
	    mean_array = numpy.zeros(dimensions[0]*dimensions[1]*dimensions[2])
	    std_array  = numpy.zeros(dimensions[0]*dimensions[1]*dimensions[2])
	    for affineFEMTranform in TransformList:
	      # get homogenius 4x4 matrix  of the form
	      #               A | b
	      #    matrix =   -----
	      #               0 | 1
	      #   
	      matrix = affineFEMTranform.GetConcatenatedTransform(0).GetMatrix()
	      #print matrix 
	      RotationMatrix = [[matrix.GetElement(0,0),matrix.GetElement(0,1),matrix.GetElement(0,2)],
				[matrix.GetElement(1,0),matrix.GetElement(1,1),matrix.GetElement(1,2)],
				[matrix.GetElement(2,0),matrix.GetElement(2,1),matrix.GetElement(2,2)]]
	      Translation =     [matrix.GetElement(0,3),matrix.GetElement(1,3),matrix.GetElement(2,3)] 
	      #print RotationMatrix 
	      print Translation 
	      TransformedFEMMesh = None
	      if vtkReflectY.GetOutput().IsA("vtkMultiBlockDataSet"):
		AppendBlocks = vtk.vtkAppendFilter()
		iter = vtkReflectY.GetOutput().NewIterator()
		iter.UnRegister(None)
		iter.InitTraversal()
		# loop over blocks...
		while not iter.IsDoneWithTraversal():
		    curInput = iter.GetCurrentDataObject()
		    vtkTransformFEMMesh = vtk.vtkTransformFilter()
		    vtkTransformFEMMesh.SetTransform( affineFEMTranform )
		    vtkTransformFEMMesh.SetInput( curInput )
		    vtkTransformFEMMesh.Update()
		    AppendBlocks.AddInput( vtkTransformFEMMesh.GetOutput() ) 
		    AppendBlocks.Update( ) 
		    iter.GoToNextItem();
		TransformedFEMMesh = AppendBlocks.GetOutput() 
	      else: 
		vtkTransformFEMMesh = vtk.vtkTransformFilter()
		vtkTransformFEMMesh.SetTransform( affineFEMTranform )
		vtkTransformFEMMesh.SetInput( vtkReflectY.GetOutput() )
		vtkTransformFEMMesh.Update()
		TransformedFEMMesh = vtkTransformFEMMesh.GetOutput() 

	      # reuse ShiftScale Geometry
	      vtkResample = vtk.vtkCompositeDataProbeFilter()
	      vtkResample.SetInput( templateImage )
	      vtkResample.SetSource( TransformedFEMMesh ) 
	      vtkResample.Update()
	      fem_point_data= vtkResample.GetOutput().GetPointData() 
	      #meantmp =  vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('Vard0Mean')) 
	      #print meantmp.max()
	      #mean_array = mean_array + vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('Vard0Mean')) 
	      #std_array  = std_array  + vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('Vard0Kurt')) 
	      #print fem_array 
	      #print type(fem_array )

	    # average
	    #mean_array = mean_array/naverage
	    #print mean_array.max()
	    #std_array  = std_array /naverage 

	    # write numpy to disk in matlab
	    #scipyio.savemat("%s/modelstats.navg%d.%04d.mat" % (work_dir,naverage,timeID), {'spacing':spacing, 'origin':origin,'Vard0Mean':mean_array,'Vard0Kurt':std_array })

	  
	    # write output
	    print "writing ", timeID, ii, work_dir
	    vtkStatsWriter = vtk.vtkDataSetWriter()
	    vtkStatsWriter.SetFileTypeToBinary()
	    vtkStatsWriter.SetFileName("%s/modelstats.navg%d.%04d.vtk" % ( work_dir,naverage,ii ))
	    vtkStatsWriter.SetInput(vtkResample.GetOutput())
	    vtkStatsWriter.Update()

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--work_dir",
                  action="store", dest="work_dir", default=None,
                  help="project data in exodus FILE  DIR/fem_stats.e to imaging", metavar = "DIR")
parser.add_option( "--data_set",
                  action="store", dest="data_set", default=None,
                  help="select data set to run ", metavar = "Int")
(options, args) = parser.parse_args()
if (options.work_dir):
  ProjectImagingMesh(options.work_dir,options.data_set)
else:
  parser.print_help()
  print options
