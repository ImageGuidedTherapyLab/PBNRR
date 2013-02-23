#==========================================================================
#
#   Copyright Insight Software Consortium
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#          http://www.apache.org/licenses/LICENSE-2.0.txt
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#==========================================================================*/

#
#  Example on the use of Dicom Reader 
#

import itk
import sys

if len(sys.argv) < 1:
    print('Usage: ' + sys.argv[0] + ' DicomDirectory')
    sys.exit(1)

#
# Reads a 3D image in with signed short (16bits/pixel) pixel type
# and save it
#
ImageType  = itk.Image.SS3

nameGenerator = itk.GDCMSeriesFileNames.New()
nameGenerator.SetUseSeriesDetails( True ) 
nameGenerator.AddSeriesRestriction("0008|0021") 

nameGenerator.SetDirectory( sys.argv[1]) 
seriesUID = nameGenerator.GetSeriesUIDs() 

for uid in seriesUID:
   # get file names
   fileNames = nameGenerator.GetFileNames( uid ) 
   print "reading:", uid
   # read
   reader = itk.ImageSeriesReader[ImageType].New()
   dicomIO = itk.GDCMImageIO.New()
   reader.SetImageIO( dicomIO )
   reader.SetFileNames( fileNames )
   reader.Update( )
   # get dictionary info
   dictionary = dicomIO.GetMetaDataDictionary()
   PrintAllKeysInDictionary = False
   if(PrintAllKeysInDictionary): 
     for key in dictionary.GetKeys():
       print key, dictionary[key]
   print dictionary['0008|103e']
   # write
   writer = itk.ImageFileWriter[ImageType].New()
   writer.SetInput( reader.GetOutput() )
   writer.SetFileName( "%s.mha" % uid );
   writer.Update() 
