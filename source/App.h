/**
  \file App.h

 */
#ifndef App_h
#define App_h

#include <G3D/G3DAll.h>

#ifdef G3D_WINDOWS
 // Compile RtAudio with windows direct sound support
#   define __WINDOWS_DS__
#endif
#include "RtAudio.h"
#include "chuck_fft.h"
#include "EWMAFrequency.h"
#include "SonicSculpturePiece.h"
#include <mutex>
/** Global array where we dump raw audio data from RtAudio */
Array<float> g_currentAudioBuffer;
int g_sampleWindowIndex;


class App : public GApp {
protected:
    G3D_DECLARE_ENUM_CLASS(AppMode, DEFAULT, MAKING_SCULPTURE);
    AppMode m_appMode;

    RealTime m_initialTime;

    
    RealTime m_lastInterestingEventTime;

    RtAudio m_rtAudio;

    struct PlayPulse {
        int initialSample;
        shared_ptr<SonicSculpturePiece> piece;
    };
    Array<PlayPulse> m_currentPlayPulses;


    struct PlayPlane {
        Vector3 direction;
        Point3 origin;
        int endWindowIndex;
        int beginWindowIndex;
      CFrame frame;
    };
    
    Array<PlayPlane> m_playPlanes;

    /** EWMA of RMS */
    float m_smoothedRootMeanSquare;

    /** 3 different rates of exponentially-weighted moving averages of frequencies */
    EWMAFrequency m_fastMovingAverage;
    EWMAFrequency m_slowMovingAverage;
    EWMAFrequency m_glacialMovingAverage;


    /** Settings for RtAudio. We never need to change the defaults */
    struct AudioSettings {
      int numChannels;
      int sampleRate;
      RtAudioFormat rtAudioFormat;
      
      AudioSettings() :
          numChannels(1),
          sampleRate(48000),
          rtAudioFormat(RTAUDIO_FLOAT32) {}
    } m_audioSettings;

    /** How many slices of time to save */
    int m_maxSavedTimeSlices;


	int m_lastSampleWindowProcessed;

    /** CPU storage of fft samples */
    Array<complex> m_cpuFrequencyAudioData;

    /** GPU storage of raw samples */
    shared_ptr<Texture> m_rawAudioTexture;
    /** GPU storage of fft samples */
    shared_ptr<Texture> m_frequencyAudioTexture;

    Array<shared_ptr<SonicSculpturePiece> > m_sonicSculpturePieces;
    shared_ptr<SonicSculpturePiece> m_currentSonicSculpturePiece;

    /** Set all our audio textures on \param Args */
    void setAudioShaderArgs(Args& args);

    /** Do all of the interaction with RtAudio that we need to to set up realtime audio capture */
    void initializeAudio();

    /** Called from onInit */
    void makeGUI();

    void handlePlayPulses();

    void generatePlayPulse(shared_ptr<SonicSculpturePiece> piece);


	void playSculpture(const Ray& playRay);

public:

    App(const GApp::Settings& settings = GApp::Settings());

    virtual void onInit() override;
    virtual void onGraphics3D(RenderDevice* rd, Array< shared_ptr<Surface> >& surface) override;
    virtual void onGraphics2D(RenderDevice* rd, Array< shared_ptr<Surface2D> >& surface2D) override;

    virtual bool onEvent(const GEvent& e) override;
    virtual void onCleanup() override;

    virtual void onUserInput(UserInput* ui) override;

    virtual void onSimulation(RealTime rdt, SimTime sdt, SimTime idt) override;

	/** Called after every audio update to get the latest audio data and compute statistics such as RMS */
	void updateAudioData();

	void updateSonicSculpture(int audioSampleOffset, int audioSampleCount);

	// TODO: rethink access patterns to not have public member variables
	/** Make sure raw audio is only touched by one thread at a time*/
	std::mutex m_rawAudioMutex;
	/** CPU storage of raw samples */
	Array<float> m_cpuRawAudioData;

};

#endif
