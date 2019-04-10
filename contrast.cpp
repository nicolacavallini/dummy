/** <pre>------------------------------------------------------------------------------------------
*                                      PPT Vision, Inc.
*                                     Copyright (c) 2001
*                                Proprietary and Confidential
*                                    All Rights Reserved
*
* $Workfile: ContrastTool1.cpp $
*
* Description:
*
* $Revision: 4 $
* $Date: 31/05/13 14:09 $
*
*                                                                                           </pre>
*/ //---------------------------------------------------------------------------------------------
#include <DisableWarning.h>
#include <ContrastTool1.h>
#ifndef toolclassinfohelper_h
#include <ToolClassInfoHelper.h>
#endif
#ifndef errortype_h
#include <ErrorType.h>
#endif
#ifndef roi_h
#include <ROI.h>
#endif
#ifndef histgram_h
#include <Histgram.h>
#endif
#ifndef hpclock_h
#include <HPClock.h>
#endif
#ifndef dbgnew_h
#include <DbgNew.h>
#endif
#ifndef camera_h
#include <Camera.h>
#endif
#ifndef _LIMITS_
#include <limits>
#endif
#ifndef vdmimage_h
#include <VdmImage.h>
#endif

ToolClassInfoHelper<ContrastTool1> ContrastTool1ClassInfo
(kContrastTool1, "Contrast", "Tool for measuring contrast in a defined area.",kFullLicense,kToolSupported,1,true);

ContrastTool1::ContrastTool1()
: ToolHelper(),
portImageIn_("Input Image"),
portShapeList_("Shape List", false, false, true),
portToolOrigin_("Tool Origin", MOrigin<double>()),
portThresholdType_("Grey Level Threshold Type"),
portFixedRange_("Fixed Threshold Range"),
portEnablePercentInRangeTest_("Enable In Range Percent Test", false),
portPercentInRangeTestRange_("Set In Range Percent"),
portEnableAreaInRangeTest_("Enable In Range Area Test", false),
portAreaInRangeTestRange_("Set In Range Area"),
portEnablePercentOutRangeTest_("Enable Out of Range Percent Test", false),
portPercentOutRangeTestRange_("Set Out of Range Percent"),
portEnableAreaOutRangeTest_("Enable Out of Range Area Test", false),
portAreaOutRangeTestRange_("Set Out of Range Area"),
portPassed_("Passed", false),
portPercentInRange_("Percent In Range", RealVal(0.0, 0.0, 100.0)), // 0~100
portPercentOutRange_("Percent Out of Range", RealVal(100.0, 0.0, 100.0)), // 0~100
portInRangeArea_("In Range Area", RealVal(0., 0., std::numeric_limits<double>::infinity())),
portOutRangeArea_("Out of Range Area", RealVal(0., 0., std::numeric_limits<double>::infinity())),
portThresholdsUsed_("Threshold Range Used"),
portHistogramSamplingStep_("Histogram Sampling Step")
{
	// register for camera callback
	Callback0 cb = MakeCallback(static_cast<Callback0*>(0), *this, &ContrastTool1::CalibrationChangedCB_);
	Camera::GetCamera()->AddCalibrationCallback(cb);
	Camera::GetCamera()->AddCameraConnectedCallback(cb);

	// add port for shape list region, and set up callbacks
	ShapeListPortDesc* pPortDesc = dynamic_cast<ShapeListPortDesc*>(portShapeList_.GetPortDesc());
	pPortDesc->SetCallbacks(MakeCallback(static_cast<Callback1wRet<ErrorType&, MPtrList<MShape*>*>*>(0), *this, &ContrastTool1::GetShapeListPatchCB_),
		MakeCallback(static_cast<Callback1wRet<MPtrList<MShape*>* const &, ErrorType>*>(0), *this, &ContrastTool1::SetShapeListPatchCB_));

	inPortList_.push_back(&portImageIn_);
	inPortList_.push_back(&portToolOrigin_);
	inPortList_.push_back(pPortDesc->AsInPort()); // shape list port
	inPortList_.push_back(&portThresholdType_);
	inPortList_.push_back(&portFixedRange_);
	inPortList_.push_back(&portEnablePercentInRangeTest_);
	inPortList_.push_back(&portPercentInRangeTestRange_);
	inPortList_.push_back(&portEnableAreaInRangeTest_);
	inPortList_.push_back(&portAreaInRangeTestRange_);
	inPortList_.push_back(&portEnablePercentOutRangeTest_);
	inPortList_.push_back(&portPercentOutRangeTestRange_);
	inPortList_.push_back(&portEnableAreaOutRangeTest_);
	inPortList_.push_back(&portAreaOutRangeTestRange_);
	inPortList_.push_back(&portHistogramSamplingStep_);
	portThresholdType_.WritePort(kAdaptiveSingleThresholdDark);
	portFixedRange_.WritePort(M1DRange<double>(0., 50.)); // default [0,50]
	portPercentInRangeTestRange_.WritePort(M1DRange<double>(0., 100.));
	portAreaInRangeTestRange_.WritePort(M1DRange<double>(0., std::numeric_limits<double>::infinity()));
	portPercentOutRangeTestRange_.WritePort(M1DRange<double>(0., 100.));
	portAreaOutRangeTestRange_.WritePort(M1DRange<double>(0., std::numeric_limits<double>::infinity()));

	outPortList_.push_back(&portPassed_);
	outPortList_.push_back(&portPercentInRange_);
	outPortList_.push_back(&portPercentOutRange_);
	outPortList_.push_back(&portInRangeArea_);
	outPortList_.push_back(&portOutRangeArea_);
	outPortList_.push_back(&portThresholdsUsed_);
	portThresholdsUsed_.SetData(M1DRange<double>(0., 50.));

	//calibratedInputPortList_.push_back(&portToolOrigin_);
	calibratedInputPortList_.push_back(pPortDesc->AsInPort());
}

ContrastTool1::~ContrastTool1()
{
	// unregister for camera callback
	Callback0 cb = MakeCallback(static_cast<Callback0*>(0), *this, &ContrastTool1::CalibrationChangedCB_);
	Camera::GetCamera()->RemoveCalibrationCallback(cb);
	Camera::GetCamera()->RemoveCameraConnectedCallback(cb);
}

void ContrastTool1::CalibrationChangedCB_()
{
	CalibrationChangedCB(portToolOrigin_.ReadPort(), roiHelper_);
}

MPtrList<MShape*>* ContrastTool1::GetShapeListPatchCB_(ErrorType& err)
{
	err = kOk;
	return(roiHelper_.GetShapeList());
}

ErrorType ContrastTool1::SetShapeListPatchCB_(MPtrList<MShape*>* const & pShapeList)
{
	SetShapeListPatchCB(pShapeList, portImageIn_.ReadPort(), portToolOrigin_.ReadPort(), roiHelper_, true);
	return(kOk);
}


ErrorType ContrastTool1::Run
(const bool& forceAbort)
{

#if defined(IMPACT_OS_DSPBIOS) && defined(WATCHDOG_ENABLED)
	MiscWatchDogRefresh();
#endif

	ErrorType eErr = kFailed;
	portPassed_.SetData(false);
	passed_ = false;

	MBaseImage* pImage = NULL;

	try
	{
		pImage = portImageIn_.ReadPort();

		if (!pImage)
			throw(kNoImage);

		ImagePixelType pixelType = pImage->GetPixelType();

		if (pixelType == kRGBTriplePixel)
			throw(kColortoGrayscaleImageConversionNeeded);

		//Getting grayscale upper bound:
		double grayUpperBound = pixelType == k8BitUIntPixel ?
			(std::numeric_limits<unsigned char>::max)()
			:
			(std::numeric_limits<unsigned short>::max)();

		// initialize contrast
		double contrast = 0.0;
		unsigned long numPixels = 0;

		// set work image
		MBaseImage* pWorkIm = 0;
		ShapeListRoiHelperHolder roiHelperHolder(roiHelper_);
		MPtrList<MBaseROI*>* pRoiList = roiHelper_.GetROI(pImage,
			portToolOrigin_.ReadPort(), eErr);


		if (eErr == kOk)
		{
			pWorkIm = pRoiList->front();
		}
		else if (eErr == kNullDataRef) // shape list pointer is null, use full image
		{
			pWorkIm = pImage;
		}
		else
			throw(eErr);

		// set threshold range
		M1DRange<double> threshRange(portFixedRange_.ReadPort());
		GreyLevelThresholdType1 threshType = static_cast<GreyLevelThresholdType1>(portThresholdType_.ReadPort().GetVal());
		switch (threshType)
		{
		case kAdaptiveSingleThresholdBright1:
		case kAdaptiveSingleThresholdDark1:
		case kAdaptiveDoubleThreshold1:
		case kAdaptiveAverageFixedWidth1:
		{
											bool isPixel16BitsDeep = pImage->GetPixelType() == k16BitUIntPixel;

											int subsamplingStep = isPixel16BitsDeep && threshType == kAdaptiveDoubleThreshold1 ?
												pow(2.0, portHistogramSamplingStep_.ReadPort())
												:
												1;

											int histogramSize = (grayUpperBound + 1.0) / subsamplingStep;

											MHistogram<double> hist(histogramSize, M1DRange<double>(0, grayUpperBound));
											pWorkIm->GetHistogram(hist, forceAbort);
											switch (threshType)
											{
											case kAdaptiveSingleThresholdBright1:
												// note: start not inclusive for bright
#ifndef IMPACT_OS_DSPBIOS
												threshRange.SetRange(hist.GetOtsuThreshold(), grayUpperBound, false, true);
#else
												threshRange.SetRange(hist.GetOtsuThreshold_Fast(), grayUpperBound, false, true);
#endif
												break;
											case kAdaptiveSingleThresholdDark1:
#ifndef IMPACT_OS_DSPBIOS
												threshRange.SetRange(0, hist.GetOtsuThreshold(), true, true);
#else
												threshRange.SetRange(0, hist.GetOtsuThreshold_Fast(), true, true);
#endif
												break;
											case kAdaptiveDoubleThreshold1:
											{
												if (isPixel16BitsDeep)
												{
													threshRange = hist.GetOtsuThresholdRangeSparse();
												}
												else
												{
													threshRange = hist.GetOtsuThresholdRange();
												}
											}
												break;
											case kAdaptiveAverageFixedWidth1:
											{
												double average = hist.GetAverage();
												double offset = 0.5*(grayUpperBound / 100.0)*threshRange.GetSize();

												threshRange.SetRange((std::max)(0., average - offset), (std::min)(grayUpperBound, average + offset));

											}
												break;
											}

											unsigned long numInRange = hist.GetAmountOfPixelsInRange(threshRange, grayUpperBound);
											numPixels = hist.GetTotalAmountOfPixels();

											if (numPixels > 0)
											    contrast = 100.0 * numInRange / numPixels;
											else
											  contrast = -1.0;


		}
			break;
		case kFixedDoubleThreshold1:
		{
									   threshRange.SetRange((std::max)(0., threshRange.GetStart()*(grayUpperBound / 100.0)), (std::min)(grayUpperBound, threshRange.GetEnd()*(grayUpperBound / 100.0)));

									   contrast = pWorkIm->GetContrast(threshRange, &numPixels, forceAbort);
		}
			break;
		default:
			throw(kInvalidParm);
		}
		portThresholdsUsed_.SetData(M1DRange<double>(threshRange.GetStart() / (grayUpperBound / 100.0), threshRange.GetEnd() / (grayUpperBound / 100.0)));

		if (AbortCheck(forceAbort))
			throw(kAborted);
		if (contrast >= 0)
		{
			portPercentInRange_.SetData(contrast);
			portPercentOutRange_.SetData(100.0 - contrast);
			portInRangeArea_.SetData(0.01 * contrast * numPixels *
				pow((pImage->GetCalibration()).GetUnitsPerPixel(), 2.0));
			portOutRangeArea_.SetData(0.01 * (100. - contrast) * numPixels *
				pow((pImage->GetCalibration()).GetUnitsPerPixel(), 2.0));
			eErr = kOk;

			bool passed = true;
			if (portEnablePercentInRangeTest_.ReadPort() && !portPercentInRangeTestRange_.ReadPort().InRange(portPercentInRange_.ReadPort()))
				passed = false;
			if (portEnableAreaInRangeTest_.ReadPort() && !portAreaInRangeTestRange_.ReadPort().InRange(portInRangeArea_.ReadPort()))
				passed = false;
			if (portEnablePercentOutRangeTest_.ReadPort() && !portPercentOutRangeTestRange_.ReadPort().InRange(portPercentOutRange_.ReadPort()))
				passed = false;
			if (portEnableAreaOutRangeTest_.ReadPort() && !portAreaOutRangeTestRange_.ReadPort().InRange(portOutRangeArea_.ReadPort()))
				passed = false;
			portPassed_.SetData(passed);
			passed_ = true;
		}

	}
	catch (MOutOfRangeError)
	{
		eErr = kRoiOffImage;
	}
	catch (std::bad_alloc &)
	{
		eErr = kOutOfMemory;
	}
	catch (ErrorType err)
	{
		eErr = err;
	}
	catch (...)
	{
		eErr = kUnhandledException;
	}

	// log error
	if (eErr != kOk)
		LogError("%s: %s", GetName().c_str(), VisionDeviceManager::GetErrorTypeString(eErr));

	return(eErr);
}

bool ContrastTool1::StreamIn
(MStream& in)
{
	// Stream in version number of this class.
	int version = 0;
	bool ok = in.Read("ver", version);
	// Switch on version number...
	switch (version)
	{
	case 1:
	case 2:
	case 3:
	{
			  // Stream in base class.
			  ok = ToolHelper::StreamIn(in) && ok;
			  // Stream in members of this class.
			  MOrigin<double> origin;
			  ok = in.Read("portToolOrigin", origin) && ok;
			  ok = (portToolOrigin_.WritePort(origin) == kOk) && ok;

			  IntVal intVal;
			  ok = in.Read("portThresholdType", intVal) && ok;
			  ok = (portThresholdType_.WritePort(intVal) == kOk) && ok;

			  M1DRange<double> range;
			  ok = in.Read("portFixedRange", range) && ok;
			  ok = (portFixedRange_.WritePort(range) == kOk) && ok;

			  MPtrList<MShape*>* pShapeList = 0;
			  ok = in.Read("shapeList", pShapeList) && ok;
			  // create shape list roi and hold on to the shape list pointer
			  SetShapeListPatchCB_(pShapeList);
			  if (pShapeList)
				  pShapeList->Release();

			  bool bVal;
			  if (version>1)
			  {
				  ok = in.Read("portEnablePercentInRangeTest", bVal) && ok;
				  ok = (portEnablePercentInRangeTest_.WritePort(bVal) == kOk) && ok;
				  ok = in.Read("portPercentInRangeTestRange", range) && ok;
				  ok = (portPercentInRangeTestRange_.WritePort(range) == kOk) && ok;
				  ok = in.Read("portEnableAreaInRangeTest", bVal) && ok;
				  ok = (portEnableAreaInRangeTest_.WritePort(bVal) == kOk) && ok;
				  ok = in.Read("portAreaInRangeTestRange", range) && ok;
				  ok = (portAreaInRangeTestRange_.WritePort(range) == kOk) && ok;
				  ok = in.Read("portEnablePercentOutRangeTest", bVal) && ok;
				  ok = (portEnablePercentOutRangeTest_.WritePort(bVal) == kOk) && ok;
				  ok = in.Read("portPercentOutRangeTestRange", range) && ok;
				  ok = (portPercentOutRangeTestRange_.WritePort(range) == kOk) && ok;
				  ok = in.Read("portEnableAreaOutRangeTest", bVal) && ok;
				  ok = (portEnableAreaOutRangeTest_.WritePort(bVal) == kOk) && ok;
				  ok = in.Read("portAreaOutRangeTestRange", range) && ok;
				  ok = (portAreaOutRangeTestRange_.WritePort(range) == kOk) && ok;
			  }
			  if (version > 2)
			  {
				  ok = in.Read("portHistogramSamplingStep", intVal) && ok;
				  ok = (portHistogramSamplingStep_.WritePort(intVal) == kOk) && ok;
			  }

			  break;
	}
	default:
	{
			   ok = false;
			   throw(MStreamError("ContrastTool1::StreamIn - bad version found!"));
			   break;
	}
	}
	return(ok);
}

bool ContrastTool1::StreamOut
(MStream& out)
const
{
	// Stream out version number of this class.
	bool ok = out.Write("ver", (int)3);

	// Stream out base class.
	ok = ToolHelper::StreamOut(out) && ok;
	// Stream out members of this class.
	ok = out.Write("portToolOrigin", portToolOrigin_.ReadPort()) && ok;
	ok = out.Write("portThresholdType", portThresholdType_.ReadPort()) && ok;
	ok = out.Write("portFixedRange", portFixedRange_.ReadPort()) && ok;
	ok = out.Write("shapeList", roiHelper_.GetShapeList()) && ok;
	ok = out.Write("portEnablePercentInRangeTest", portEnablePercentInRangeTest_.ReadPort()) && ok;
	ok = out.Write("portPercentInRangeTestRange", portPercentInRangeTestRange_.ReadPort()) && ok;
	ok = out.Write("portEnableAreaInRangeTest", portEnableAreaInRangeTest_.ReadPort()) && ok;
	ok = out.Write("portAreaInRangeTestRange", portAreaInRangeTestRange_.ReadPort()) && ok;
	ok = out.Write("portEnablePercentOutRangeTest", portEnablePercentOutRangeTest_.ReadPort()) && ok;
	ok = out.Write("portPercentOutRangeTestRange", portPercentOutRangeTestRange_.ReadPort()) && ok;
	ok = out.Write("portEnableAreaOutRangeTest", portEnableAreaOutRangeTest_.ReadPort()) && ok;
	ok = out.Write("portAreaOutRangeTestRange", portAreaOutRangeTestRange_.ReadPort()) && ok;
	ok = out.Write("portHistogramSamplingStep", portHistogramSamplingStep_.ReadPort()) && ok;

	return(ok);
}
