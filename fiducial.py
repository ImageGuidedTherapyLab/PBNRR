
RedViewer = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeRed')
#RedViewer.SetOrientationToCoronal()
RedViewer.SetOrientationToAxial()


fidID = "vtkMRMLAnnotationFiducialNode1"
fid = slicer.mrmlScene.GetNodeByID(fidID)
coords = [0,0,0]
fid.GetFiducialCoordinates(coords)
print coords[0], coords[1], coords[2]

quaterion = RedViewer.GetSliceToRAS()
print quaterion
quaterion.SetElement(2,3,coords[2])
RedViewer.SetSliceToRAS(quaterion)
RedViewer.Modified()
slicer.app.applicationLogic().PropagateVolumeSelection(0)
print RedViewer.GetSliceToRAS()
##     ijkToRAS = vtk.vtkMatrix4x4()
##     inputVolume.GetIJKToRASMatrix(ijkToRAS)
##     outputVolume.SetIJKToRASMatrix(ijkToRAS)
