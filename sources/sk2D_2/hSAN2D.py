try: paraview.simple
except: from paraview.simple import *

import sys

i = int(sys.argv[1])

paraview.simple._DisableFirstRenderCameraReset()

RenderView1 = GetRenderView()
my_2d1_vts = XMLStructuredGridReader( FileName=['my_2d%i.vts' % (i)] )

RenderView1.CenterAxesVisibility = 0
RenderView1.OrientationAxesVisibility = 0

my_2d1_vts.PointArrayStatus = ['Unnamed0']

DataRepresentation1 = Show()
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]

RenderView1.CenterOfRotation = [65.5, 82.5, 0.0]

Transform1 = Transform( Transform="Transform" )

a1_Unnamed0_PVLookupTable = GetLookupTableForArray( "Unnamed0", 1, NanColor=[0.25, 0.0, 0.0], RGBPoints=[0.0, 0.23, 0.299, 0.754, 1.0, 0.7059953924945706, 0.01625293517428646, 0.1500022472819457], VectorMode='Magnitude', ColorSpace='Diverging', LockScalarRange=1 )

a1_Unnamed0_PiecewiseFunction = CreatePiecewiseFunction( Points=[-0.1, 0.0, 1.0, 1.0] )

RenderView1.CameraPosition = [65.5, 82.5, 407.0022200734933]
RenderView1.CameraFocalPoint = [65.5, 82.5, 0.0]
RenderView1.CameraClippingRange = [402.93219787275837, 413.10725337459564]
RenderView1.CameraParallelScale = 105.33992595402752

DataRepresentation1.ColorArrayName = 'Unnamed0'
DataRepresentation1.LookupTable = a1_Unnamed0_PVLookupTable

Transform1.Transform = "Transform"

DataRepresentation2 = Show()
DataRepresentation2.ColorArrayName = 'Unnamed0'
DataRepresentation2.LookupTable = a1_Unnamed0_PVLookupTable
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]

DataRepresentation1.Visibility = 0

Transform1.Transform.Scale = [1.0, -1.0, 0.0]

RenderView1.CameraViewUp = [0.0, 1.0, 0.0]
RenderView1.CameraPosition = [65.5, -82.5, 407.0022200734933]
RenderView1.CameraFocalPoint = [65.5, -82.5, 0.0]
RenderView1.CenterOfRotation = [65.5, -82.5, 0.0]
RenderView1.CameraClippingRange = [402.93219787275837, 413.10725337459564]
WriteImage('reentry2d%06i.png' % (i))
Render()


