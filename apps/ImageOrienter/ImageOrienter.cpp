/*=========================================================================

Program: Image reorienter
Language: C++
Date: $Date: 2020-09-04 $
Version: $Revision: 1 $
Author: $Ricardo Corredor$

==========================================================================*/

#include "itkImage.h"
#include "itkSpatialOrientationAdapter.h"
#include "itkOrientImageFilter.h"
#include "itkChangeInformationImageFilter.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"

#include "itkParameterFileParser.h"

#define ORIENTATION_ASL itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL


#include <metaCommand.h>


template<typename TImagePixelType>
void ProcessImage(std::string p_inputImgFilePath, std::string p_outputImgFilePath, std::string p_imgParamsFilePath, int p_mode, int verbose)
{
	const unsigned int dimension = 3;
	typedef itk::Image<TImagePixelType, dimension> InputImageType;
	typedef itk::Image<TImagePixelType, dimension> OutputImageType;
	typedef itk::ImageFileReader<InputImageType> ReaderType;
	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	typedef itk::SpatialOrientationAdapter SpatialOrientationAdapterType;
	typedef itk::OrientImageFilter<InputImageType, InputImageType> OrientInputImageFilterType;

	if (verbose != 0)
		std::cerr << "Read input image... " << std::endl;

	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(p_inputImgFilePath);
	reader->Update();

	typename InputImageType::Pointer inputImage = reader->GetOutput();

	typename InputImageType::PointType inputImgOrigin = inputImage->GetOrigin();
	typename InputImageType::SpacingType inputImgSpacing = inputImage->GetSpacing();
	typename InputImageType::DirectionType inputImgDirection = inputImage->GetDirection();
	SpatialOrientationAdapterType orientationAdapter;
	SpatialOrientationAdapterType::OrientationType inputImgOrientation = orientationAdapter.FromDirectionCosines(inputImgDirection);

	if (verbose != 0)
	{
		std::cout << "Original origin: ";
		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << inputImgOrigin[0] << " " << inputImgOrigin[1] << " " << inputImgOrigin[2] << std::endl;
		std::cout << "Original spacing: ";
		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << inputImgSpacing[0] << " " << inputImgSpacing[1] << " " << inputImgSpacing[2] << std::endl;
		std::cout << "Original direction:" << std::endl;
		std::cout << inputImgDirection << std::endl;
		std::cout << "Original orientation: ";
		std::cout << inputImgOrientation << std::endl;
	}

	if (p_mode == 0) //Prepare for preprocessing
	{
		typename OrientInputImageFilterType::Pointer orienter = OrientInputImageFilterType::New();
		orienter->SetGivenCoordinateOrientation(inputImgOrientation);
		orienter->SetInput(inputImage);
		orienter->SetDesiredCoordinateOrientation(ORIENTATION_ASL);
		orienter->Update();
		typename InputImageType::Pointer inputImageASL = orienter->GetOutput();

		typename InputImageType::PointType newImgOrigin = inputImageASL->GetOrigin();

		newImgOrigin[0] = 0;
		newImgOrigin[1] = 0;
		newImgOrigin[2] = 0;


		typename InputImageType::SpacingType newImgSpacing = inputImageASL->GetSpacing();

		//  Set the spacing of the image
		newImgSpacing[0] = floor(newImgSpacing[0] * 100.0 + 0.5) / 100.0;
		newImgSpacing[1] = floor(newImgSpacing[1] * 100.0 + 0.5) / 100.0;
		newImgSpacing[2] = floor(newImgSpacing[2] * 100.0 + 0.5) / 100.0;


		typename InputImageType::DirectionType newImgDirection = inputImageASL->GetDirection();
		newImgDirection.SetIdentity();
		newImgDirection[0][0] = 0;
		newImgDirection[1][1] = 0;
		newImgDirection[2][2] = 0;
		newImgDirection[0][2] = -1;
		newImgDirection[1][0] = 1;
		newImgDirection[2][1] = -1;

		if (verbose != 0)
		{
			std::cout << "New origin: ";
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << newImgOrigin[0] << " " << newImgOrigin[1] << " " << newImgOrigin[2] << std::endl;
			std::cout << "New spacing: ";
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << newImgSpacing[0] << " " << newImgSpacing[1] << " " << newImgSpacing[2] << std::endl;
			std::cout << "New direction:" << std::endl;
			std::cout << newImgDirection << std::endl;
		}

		typedef itk::ChangeInformationImageFilter<InputImageType> RoundSpcFilterType;
		typename RoundSpcFilterType::Pointer changeImgPropsFilter = RoundSpcFilterType::New();
		changeImgPropsFilter->SetOutputSpacing(newImgSpacing);
		changeImgPropsFilter->SetOutputDirection(newImgDirection);
		changeImgPropsFilter->SetOutputOrigin(newImgOrigin);
		changeImgPropsFilter->SetChangeOrigin(true);
		changeImgPropsFilter->SetChangeSpacing(true);
		changeImgPropsFilter->SetChangeDirection(true);

		try {
			changeImgPropsFilter->SetInput(inputImageASL);
			changeImgPropsFilter->Update();
		}
		catch (itk::ExceptionObject exc)
		{
			std::cerr << "Error:" << std::endl;
			std::cerr << exc << std::endl;
		}

		typename OutputImageType::Pointer outputImage = changeImgPropsFilter->GetOutput();

		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(p_outputImgFilePath);
		writer->SetInput(outputImage);
		writer->Update();

		std::ofstream propertiesfile;
		propertiesfile.open(p_imgParamsFilePath);
		propertiesfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "(ImageOrigin " << inputImgOrigin[0] << " " << inputImgOrigin[1] << " " << inputImgOrigin[2] << ")\n";
		propertiesfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "(ImageSpacing " << inputImgSpacing[0] << " " << inputImgSpacing[1] << " " << inputImgSpacing[2] << ")\n";

		propertiesfile << "(ImageDirection ";
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
				propertiesfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << inputImgDirection[j][i] << " ";
			}
		}
		propertiesfile << ")\n";

		propertiesfile.close();
	}
	else //Recover image properties
	{
		typedef itk::ParameterFileParser ParameterFileParserType;
		typename ParameterFileParserType::Pointer paramsParser = ParameterFileParserType::New();
		paramsParser->SetParameterFileName(p_imgParamsFilePath);
		try
		{
			paramsParser->ReadParameterFile();
			typename ParameterFileParserType::ParameterMapType propertiesMap = paramsParser->GetParameterMap();
			std::vector<std::string> imgOriginArray = propertiesMap["ImageOrigin"];
			std::vector<std::string> imgSpacingArray = propertiesMap["ImageSpacing"];
			std::vector<std::string> imgDirectionArray = propertiesMap["ImageDirection"];

			if (verbose != 0)
			{
				std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "origin[0] = " << imgOriginArray[0] << std::endl;
				std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "origin[1] = " << imgOriginArray[1] << std::endl;
				std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "origin[2] = " << imgOriginArray[2] << std::endl;

				std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "spacing[0] = " << imgSpacingArray[0] << std::endl;
				std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "spacing[1] = " << imgSpacingArray[1] << std::endl;
				std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "spacing[2] = " << imgSpacingArray[2] << std::endl;
			}

			typename InputImageType::DirectionType tmpDirection; // = tmpImg->GetDirection();
			//tmpDirection.SetIdentity();
			if (verbose != 0)
				std::cout << "direction ";
			int count = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					if (verbose != 0)
						std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << imgDirectionArray[count] << " ";

					tmpDirection[i][j] = std::stof(imgDirectionArray[count]);
					count++;
				}
			}

			if (verbose != 0)
				std::cout << "" << std::endl;

			SpatialOrientationAdapterType orientationAdapter;
			typename SpatialOrientationAdapterType::OrientationType imageOrientation = orientationAdapter.FromDirectionCosines(tmpDirection);

			if (verbose != 0)
				std::cout << "Orientation:" << imageOrientation << std::endl;

			typename OrientInputImageFilterType::Pointer orienter = OrientInputImageFilterType::New();
			orienter->SetGivenCoordinateOrientation(ORIENTATION_ASL);
			orienter->SetInput(inputImage);
			orienter->SetDesiredCoordinateOrientation(imageOrientation);
			orienter->Update();
			typename InputImageType::Pointer inputImageRecovered = orienter->GetOutput();

			typename InputImageType::PointType newImgOrigin;

			newImgOrigin[0] = std::stof(imgOriginArray[0]);
			newImgOrigin[1] = std::stof(imgOriginArray[1]);
			newImgOrigin[2] = std::stof(imgOriginArray[2]);

			typename InputImageType::SpacingType newImgSpacing;

			//  Set the spacing of the image
			newImgSpacing[0] = std::stof(imgSpacingArray[0]);
			newImgSpacing[1] = std::stof(imgSpacingArray[1]);
			newImgSpacing[2] = std::stof(imgSpacingArray[2]);


			typedef itk::ChangeInformationImageFilter<InputImageType> RoundSpcFilterType;
			typename RoundSpcFilterType::Pointer changeImgPropsFilter = RoundSpcFilterType::New();
			changeImgPropsFilter->SetOutputSpacing(newImgSpacing);
			changeImgPropsFilter->SetOutputDirection(tmpDirection);
			changeImgPropsFilter->SetOutputOrigin(newImgOrigin);
			changeImgPropsFilter->SetChangeOrigin(true);
			changeImgPropsFilter->SetChangeSpacing(true);
			changeImgPropsFilter->SetChangeDirection(true);

			try {
				changeImgPropsFilter->SetInput(inputImageRecovered);
				changeImgPropsFilter->Update();
			}
			catch (itk::ExceptionObject exc)
			{
				std::cerr << "Error:" << std::endl;
				std::cerr << exc << std::endl;
			}

			typename OutputImageType::Pointer outputImage = changeImgPropsFilter->GetOutput();

			typename WriterType::Pointer writer = WriterType::New();
			writer->SetFileName(p_outputImgFilePath);
			writer->SetInput(outputImage);
			writer->Update();

		}
		catch (itk::ExceptionObject & e)
		{
			std::cerr << "Parameters/Properties file not suppported"<< e << std::endl;

		}

	}

}



int main(int argc, char *argv[])
{

	// Process some command-line arguments intended for BatchMake
  MetaCommand command;

	command.SetOption("arg_InputFilePath", "i", true, "Absolute path of the input image");
	command.AddOptionField("arg_InputFilePath", "filename", MetaCommand::STRING, true);
	command.SetOption("arg_OutputFilePath", "o", true, "Absolute path of the output image");
	command.AddOptionField("arg_OutputFilePath", "filename", MetaCommand::STRING, true);
	command.SetOption("arg_ParametersFilePath", "p", true, "Absolute path of the image parameters file to be used");
	command.AddOptionField("arg_ParametersFilePath", "filename", MetaCommand::STRING, true);

	// Option for setting the option mode
	command.SetOption("mode", "m", true, ". 0:Prepare for preprocessing, 1:Recover image properties");
	command.SetOptionLongTag("mode", "mode");
	command.AddOptionField("mode", "value", MetaCommand::INT, true);

	// Option for setting the option verbose
	command.SetOption("verbose", "t", false, "Verbose 0 (default), 1: Print messages");
	command.SetOptionLongTag("verbose", "verbose");
	command.AddOptionField("verbose", "value", MetaCommand::INT, true);

	command.SetParseFailureOnUnrecognizedOption( true );


	if( !command.Parse( argc, argv ) )
	{
		std::cerr << "Error during " << argv[0]
		          << " command argument parsing." << std::endl;
		return EXIT_FAILURE;
	}

	std::string inputImgFilePath, imgParamsFilePath, outputImgFilePath;
	inputImgFilePath = command.GetValueAsString("arg_InputFilePath","filename");
	outputImgFilePath = command.GetValueAsString("arg_OutputFilePath","filename");
	imgParamsFilePath = command.GetValueAsString("arg_ParametersFilePath","filename");
	int mode = command.GetValueAsInt("mode", "value");

	int verbose = 0;
	if (command.GetOptionWasSet("verbose"))
		verbose = command.GetValueAsInt("verbose", "value");

	if (verbose != 0)
	{
		std::cout << "inputImgFilePath:" << inputImgFilePath << std::endl;
		std::cout << "outputImgFilePath:" << outputImgFilePath << std::endl;
		std::cout << "imgParamsFilePath:" << imgParamsFilePath << std::endl;
		std::cout << "mode:" << mode << std::endl;
	}

	typename itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputImgFilePath.c_str(), itk::CommonEnums::IOFileMode::ReadMode);
	if (!imageIO)
	{
		std::cerr << "Could not CreateImageIO for: " << inputImgFilePath << std::endl;
		return EXIT_FAILURE;
	}
	imageIO->SetFileName(inputImgFilePath);
	imageIO->ReadImageInformation();

	// itk::ImageIOBase::IOComponentEnum pixelType;
	using ScalarPixelType = itk::ImageIOBase::IOComponentEnum;
	const ScalarPixelType pixelType = imageIO->GetComponentType();
	const size_t numDimensions = imageIO->GetNumberOfDimensions();
	if (verbose != 0)
	{
		std::cout << "Pixel Type is " << imageIO->GetComponentTypeAsString(pixelType) // 'double'
			<< std::endl;
		std::cout << "numDimensions: " << numDimensions << std::endl; // '3'
	}

	switch (pixelType)
	{
	case itk::ImageIOBase::IOComponentEnum::DOUBLE:
		typedef double ImagePixelTypeDouble;
		ProcessImage<ImagePixelTypeDouble>(inputImgFilePath, outputImgFilePath, imgParamsFilePath, mode, verbose);
		break;
	case itk::ImageIOBase::IOComponentEnum::USHORT:
		typedef unsigned short ImagePixelTypeUShort;
		ProcessImage<ImagePixelTypeUShort>(inputImgFilePath, outputImgFilePath, imgParamsFilePath, mode, verbose);
		break;
	case itk::ImageIOBase::IOComponentEnum::SHORT:
		typedef unsigned short ImagePixelTypeShort;
		ProcessImage<ImagePixelTypeShort>(inputImgFilePath, outputImgFilePath, imgParamsFilePath, mode, verbose);
		break;
	case itk::ImageIOBase::IOComponentEnum::FLOAT:
		typedef float ImagePixelTypeFloat;
		ProcessImage<ImagePixelTypeFloat>(inputImgFilePath, outputImgFilePath, imgParamsFilePath, mode, verbose);
		break;
	case itk::ImageIOBase::IOComponentEnum::INT:
		typedef int ImagePixelTypeInt;
		ProcessImage<ImagePixelTypeInt>(inputImgFilePath, outputImgFilePath, imgParamsFilePath, mode, verbose);
		break;
	case itk::ImageIOBase::IOComponentEnum::UCHAR:
		typedef unsigned char ImagePixelTypeUChar;
		ProcessImage<ImagePixelTypeUChar>(inputImgFilePath, outputImgFilePath, imgParamsFilePath, mode, verbose);
		break;
	case itk::ImageIOBase::IOComponentEnum::CHAR:
		typedef unsigned char ImagePixelTypeChar;
		ProcessImage<ImagePixelTypeChar>(inputImgFilePath, outputImgFilePath, imgParamsFilePath, mode, verbose);
		break;

	default:
		std::cerr << "Pixel Type ("
			<< imageIO->GetComponentTypeAsString(pixelType)
			<< ") not supported. Exiting." << std::endl;
		return -1;
	}


  	std::cerr << "Done! " << std::endl;
  	return 0;
}
