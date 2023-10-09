/*=========================================================================

Program: N4 MRI Bias Field Correction, based on Slicer implementation
Language: C++
Date: $Date: 2019-01-28 $
Version: $Revision: 2 $

==========================================================================*/

#include "itkImage.h"
#include "itkShrinkImageFilter.h"
#include "itkExtractImageFilter.h"
#include <itkN4BiasFieldCorrectionImageFilter.h>
#include "itkOtsuThresholdImageFilter.h"

#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMultiThreaderBase.h"
#include <metaCommand.h>

int main(int argc, char *argv[])
{

	// Process some command-line arguments intended for BatchMake
  	MetaCommand command;

	command.SetOption("SOURCE_filename", "i", true, "the input image");
	command.AddOptionField("SOURCE_filename", "filename", MetaCommand::STRING, true);
	command.SetOption("MASK_filename", "m", false, " the mask");
	command.AddOptionField("MASK_filename", "filename", MetaCommand::STRING, true);
	command.SetOption("OUTPUT_filename", "o", true, "the filename of the output image");
	command.AddOptionField("OUTPUT_filename", "filename", MetaCommand::STRING, true);


	// Option for setting the ConvergenceThreshold
	command.SetOption("convergenceThreshold","c",false, "ConvergenceThreshold");
	command.SetOptionLongTag("convergenceThreshold","convergence-threshold");
	command.AddOptionField("convergenceThreshold","value",MetaCommand::FLOAT,true);

	// Option for setting the FullWidthAtHalfMaximum
	command.SetOption("fullWidthAtHalfMaximum","l",false, "FullWidthAtHalfMaximum");
	command.SetOptionLongTag("fullWidthAtHalfMaximum","full-width-half-max");
	command.AddOptionField("fullWidthAtHalfMaximum","value",MetaCommand::FLOAT,true);

	// Option for setting the WienerFilterNoise
	command.SetOption("wienerFilterNoise","w",false, "WienerFilterNoise");
	command.SetOptionLongTag("wienerFilterNoise","wiener-noise");
	command.AddOptionField("wienerFilterNoise","value",MetaCommand::FLOAT,true);

	// Option for setting the shrink factor
	command.SetOption("shrinkFactor","s",false, "Shrink factor");
	command.SetOptionLongTag("shrinkFactor","shrink-factor");
	command.AddOptionField("shrinkFactor","value",MetaCommand::FLOAT,true);

	// Option for setting the NumberOfFittingLevels
	command.SetOption("numOfFittingLevels","f",false, "NumberOfFittingLevels");
	command.SetOptionLongTag("numOfFittingLevels","num-fitting-levels");
	command.AddOptionField("numOfFittingLevels","value",MetaCommand::INT,true);

	// Option for setting the NumberOfHistogramBins
	command.SetOption("numHistogramBins","h",false, "NumberOfHistogramBins");
	command.SetOptionLongTag("numHistogramBins","num-histogram-bins");
	command.AddOptionField("numHistogramBins","value",MetaCommand::INT,true);

	// Option for setting the maximumNumberOfIterations
	command.SetOption("maxNumIterations","t",false, "maximumNumberOfIterations");
	command.SetOptionLongTag("maxNumIterations","max-num-iterations");
	command.AddOptionField("maxNumIterations","value",MetaCommand::LIST,true);	

	// Option for setting the NumThreads
	command.SetOption("numThreads", "d", false, "NumThreads");
	command.SetOptionLongTag("numThreads", "num-threads");
	command.AddOptionField("numThreads", "value", MetaCommand::INT, true);

	unsigned int ShrinkFactor = 0;
	unsigned int NumberOfFittingLevels = 0;
	unsigned int NumberOfHistogramBins = 0;
	float WienerFilterNoise = 0;
	float FullWidthAtHalfMaximum = 0;
	float ConvergenceThreshold = 0;
	int NumThreads = 0;

	std::list< std::string > maximumNumberOfIterations;
  	maximumNumberOfIterations.clear();

   	command.SetParseFailureOnUnrecognizedOption( true );


	if( !command.Parse( argc, argv ) )
	{
		std::cerr << "Error during " << argv[0]
		          << " command argument parsing." << std::endl;
		return EXIT_FAILURE;
	}

	if( command.GetOptionWasSet("convergenceThreshold") )
		ConvergenceThreshold = command.GetValueAsFloat("convergenceThreshold","value");
	else
		ConvergenceThreshold = 0.0001;

	if( command.GetOptionWasSet("fullWidthAtHalfMaximum") )
		FullWidthAtHalfMaximum = command.GetValueAsFloat("fullWidthAtHalfMaximum","value");
	else
		FullWidthAtHalfMaximum = 0.15;

	if( command.GetOptionWasSet("wienerFilterNoise") )
		WienerFilterNoise = command.GetValueAsFloat("wienerFilterNoise","value");
	else
		WienerFilterNoise = 0.01;

	if( command.GetOptionWasSet("shrinkFactor") )
		ShrinkFactor = command.GetValueAsFloat("shrinkFactor","value");
	else
		ShrinkFactor = 4;

	if( command.GetOptionWasSet("numOfFittingLevels") )
		NumberOfFittingLevels = command.GetValueAsInt("numOfFittingLevels","value");
	else
		NumberOfFittingLevels = 3;

	if( command.GetOptionWasSet("numHistogramBins") )
		NumberOfHistogramBins = command.GetValueAsInt("numHistogramBins","value");
	else
		NumberOfHistogramBins = 200;

	if( command.GetOptionWasSet("maxNumIterations") )
		maximumNumberOfIterations = command.GetValueAsList("maxNumIterations");
	else
	{
		maximumNumberOfIterations.push_back("50");
		maximumNumberOfIterations.push_back("40");
		maximumNumberOfIterations.push_back("30");
	}

	if (command.GetOptionWasSet("numThreads"))
		NumThreads = command.GetValueAsInt("numThreads", "value");
	else
		NumThreads = 4;

	itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(NumThreads);
	itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(NumThreads);

	std::string SOURCE_filename, MASK_filename, OUTPUT_filename;
	SOURCE_filename = command.GetValueAsString("SOURCE_filename","filename");	
	OUTPUT_filename = command.GetValueAsString("OUTPUT_filename","filename");

	bool isMaskProvided = false;
	if( command.GetOptionWasSet("MASK_filename") )
	{
		MASK_filename = command.GetValueAsString("MASK_filename","filename");
		isMaskProvided = true;
		std::cout << "MASK_filename:" << MASK_filename << std::endl;
 	}

	std::cout << "SOURCE_filename:" << SOURCE_filename << std::endl;
	std::cout << "OUTPUT_filename:" << OUTPUT_filename << std::endl;
	std::cout << "ShrinkFactor:" << ShrinkFactor << std::endl;
	std::cout << "NumberOfFittingLevels:" << NumberOfFittingLevels << std::endl;
	std::cout << "NumberOfHistogramBins:" << NumberOfHistogramBins << std::endl;
	std::cout << "WienerFilterNoise:" << WienerFilterNoise << std::endl;
	std::cout << "FullWidthAtHalfMaximum:" << FullWidthAtHalfMaximum << std::endl;
	std::cout << "ConvergenceThreshold:" << ConvergenceThreshold << std::endl;
	


	const unsigned int dimension = 3;   
	typedef float InputPixelType;
        typedef float OutputPixelType;
	typedef float MaskPixelType;
	typedef itk::Image<InputPixelType, dimension> InputImageType;
	typedef itk::Image<OutputPixelType, dimension> OutputImageType;
	typedef itk::Image<MaskPixelType, dimension> MaskType;    
	typedef itk::ImageFileReader<InputImageType> ReaderType;
	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	typedef itk::ImageFileReader<MaskType> MaskReaderType;


	std::cerr << "Read input image... " << std::endl;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(SOURCE_filename);
	reader->Update();
  
  	InputImageType::Pointer inputImage = reader->GetOutput();

	MaskType::Pointer maskImage = nullptr;
	if (!isMaskProvided)
	{
		std::cout << "Mask not provided.  Creaing mask based on Otsu ..." << std::endl;
		typedef itk::OtsuThresholdImageFilter<InputImageType, MaskType> ThresholderType;
		ThresholderType::Pointer otsu = ThresholderType::New();
		otsu->SetInput(inputImage);
		otsu->SetNumberOfHistogramBins(200);
		otsu->SetInsideValue(0);
		otsu->SetOutsideValue(1);
		otsu->Update();
		std::cout << "Mask created!" << std::endl;
		maskImage = otsu->GetOutput();
	}
	else
	{
		std::cout << "Read input mask:" << MASK_filename << std::endl;
		MaskReaderType::Pointer maskReader = MaskReaderType::New();
		maskReader->SetFileName(MASK_filename);
		maskReader->Update();
		maskImage = maskReader->GetOutput();
	}

	maskImage->SetOrigin(inputImage->GetOrigin());
	maskImage->SetSpacing(inputImage->GetSpacing());


	std::cerr << "Shrink input image... " << std::endl;
	typedef itk::ShrinkImageFilter<InputImageType, InputImageType> ShrinkerType;
	ShrinkerType::Pointer shrinker = ShrinkerType::New();
	  
	shrinker->SetInput(inputImage);
	shrinker->SetShrinkFactors(ShrinkFactor);
	shrinker->Update();
	shrinker->UpdateLargestPossibleRegion();


	std::cerr << "Shrink input mask... " << std::endl;
	typedef itk::ShrinkImageFilter<MaskType, MaskType> MaskShrinkerType;
	MaskShrinkerType::Pointer maskShrinker = MaskShrinkerType::New();
	  
	maskShrinker->SetInput(maskImage);
	maskShrinker->SetShrinkFactors(ShrinkFactor);
	maskShrinker->Update();
	maskShrinker->UpdateLargestPossibleRegion();

	std::cerr << "Run N4 Bias Field Correction... " << std::endl;
	typedef itk::N4BiasFieldCorrectionImageFilter<InputImageType,MaskType,InputImageType> CorrecterType;
  	CorrecterType::Pointer correcter = CorrecterType::New();  
  	correcter->SetInput1(shrinker->GetOutput()); 
  	correcter->SetMaskImage(maskShrinker->GetOutput());
  	correcter->SetMaskLabel(1);

  	//With B-spline grid res. = [1, 1, 1]
  	CorrecterType::ArrayType NumberOfControlPoints(NumberOfFittingLevels);
  	NumberOfControlPoints[0] = NumberOfFittingLevels+1;
  	NumberOfControlPoints[1] = NumberOfFittingLevels+1;
  	NumberOfControlPoints[2] = NumberOfFittingLevels+1;

	CorrecterType::VariableSizeArrayType iters(NumberOfFittingLevels); 
	std::list< std::string >::const_iterator baselineItr;
	baselineItr = maximumNumberOfIterations.begin();
	int i = 0;
	std::cout << "maximumNumberOfIterations:[";
	while (baselineItr != maximumNumberOfIterations.end())
	{
		iters[i] = std::atoi(baselineItr->c_str());
		std::cout << " " << baselineItr->c_str();
		i = i + 1;
		++baselineItr;
	}
	std::cout << "]" << std::endl;

	correcter->SetMaximumNumberOfIterations(iters);
	correcter->SetNumberOfFittingLevels(NumberOfFittingLevels);
	correcter->SetNumberOfControlPoints(NumberOfControlPoints); 
	correcter->SetWienerFilterNoise(WienerFilterNoise);
	correcter->SetBiasFieldFullWidthAtHalfMaximum(FullWidthAtHalfMaximum);
	correcter->SetConvergenceThreshold(ConvergenceThreshold);
	correcter->SetNumberOfHistogramBins(NumberOfHistogramBins);
	correcter->Update();

	std::cerr << "Extract log field estimated... " << std::endl;
	typedef CorrecterType::BiasFieldControlPointLatticeType PointType;
	typedef CorrecterType::ScalarImageType ScalarImageType;

	typedef itk::BSplineControlPointImageFilter<PointType, ScalarImageType> BSplinerType;
	BSplinerType::Pointer bspliner = BSplinerType::New();
	bspliner->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
	bspliner->SetSplineOrder( correcter->GetSplineOrder() );
	bspliner->SetSize( inputImage->GetLargestPossibleRegion().GetSize() );
	bspliner->SetOrigin( inputImage->GetOrigin() );
	bspliner->SetDirection( inputImage->GetDirection() );
	bspliner->SetSpacing( inputImage->GetSpacing() );
	bspliner->Update();

	InputImageType::Pointer logField = InputImageType::New();
	logField->SetOrigin(bspliner->GetOutput()->GetOrigin());
	logField->SetSpacing(bspliner->GetOutput()->GetSpacing());
	logField->SetRegions(bspliner->GetOutput()->GetLargestPossibleRegion().GetSize());
	logField->SetDirection(bspliner->GetOutput()->GetDirection());
	logField->Allocate();

	itk::ImageRegionIterator<ScalarImageType> ItB(bspliner->GetOutput(),bspliner->GetOutput()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<InputImageType> ItF(logField,logField->GetLargestPossibleRegion());
  
	for(ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF)
	{
		ItF.Set( ItB.Get()[0] );
	}

	std::cerr << "Compute the image corrected... " << std::endl;
	typedef itk::ExpImageFilter<InputImageType, InputImageType> ExpFilterType;
	ExpFilterType::Pointer expFilter = ExpFilterType::New();
	expFilter->SetInput( logField );
	expFilter->Update();

	typedef itk::DivideImageFilter<InputImageType, InputImageType, InputImageType> DividerType;
	DividerType::Pointer divider = DividerType::New();
	divider->SetInput1( inputImage );
	divider->SetInput2( expFilter->GetOutput() );
	divider->Update();

	InputImageType::IndexType inputImageIndex = inputImage->GetLargestPossibleRegion().GetIndex();
	InputImageType::SizeType inputImageSize = inputImage->GetLargestPossibleRegion().GetSize();

	InputImageType::RegionType inputRegion;
	inputRegion.SetIndex(inputImageIndex);
	inputRegion.SetSize(inputImageSize);

	typedef itk::ExtractImageFilter<InputImageType, InputImageType> CropperType;
	CropperType::Pointer cropper = CropperType::New();
	cropper->SetInput(divider->GetOutput());
	cropper->SetExtractionRegion(inputRegion);
	cropper->SetDirectionCollapseToSubmatrix();
	cropper->Update();

	std::cerr << "Write ouput images ... " << std::endl;

	WriterType::Pointer writer = WriterType::New();
  	writer->SetFileName(OUTPUT_filename);
  	writer->SetInput(cropper->GetOutput());
  	writer->Update();

	/*WriterType::Pointer fieldWriter = WriterType::New();
  	fieldWriter->SetFileName(argv[4]);
  	fieldWriter->SetInput(logField);
	fieldWriter->Update();*/

  	std::cerr << "Done! " << std::endl;
  	return 0;
}



