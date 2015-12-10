/** \file App.cpp */
#include "App.h"
#include "util.h"
#include "Synthesizer.h"
// Tells C++ to invoke command-line main() function even on OS X and Win32.
G3D_START_AT_MAIN();


static App* s_app = NULL;
#define DEV_MODE 1
int main(int argc, const char* argv[]) {
    (void)argc; (void)argv;
    GApp::Settings settings(argc, argv);
    
    // Change the window and other startup parameters by modifying the
    // settings class.  For example:
    settings.window.width       = 1280; 
    settings.window.height      = 720;
    settings.window.asynchronous = false;
    settings.window.resizable = true;

    settings.window.fullScreen = DEV_MODE == 0;

    settings.window.caption = "Sonic Sculpting";
    settings.dataDir = FileSystem::currentDirectory();
    settings.screenshotDirectory = "../journal/";

    

    settings.renderer.deferredShading = false;
    settings.renderer.orderIndependentTransparency = true;


    return App(settings).run();
}

// Helper function
static shared_ptr<UniversalMaterial> getSonicSculptureMaterial(int index) {
    shared_ptr<Texture> lambertianTex = Texture::createEmpty(format("Sonic Sculpture %d", index), 512, 1, ImageFormat::RGBA16F());
    static shared_ptr<Framebuffer> fb = Framebuffer::create("Sonic Sculpture Lambertian FB Clearer");
    fb->set(Framebuffer::COLOR0, lambertianTex);
    RenderDevice* rd = RenderDevice::current;
    rd->push2D(fb); {
        rd->setColorClearValue(Color3::white() * 0.9f);
        rd->clear();
    } rd->pop2D();
    lambertianTex->generateMipMaps();
    UniversalMaterial::Specification spec;
    spec.setLambertian(lambertianTex);
    static uint32 dummyBytes[512];
    for (int i = 0; i < 512; ++i) {
        dummyBytes[i] = 4294967295;
    }
    shared_ptr<Texture> emissiveTex = Texture::fromMemory(format("Sonic Sculpture %d Emissive", index), dummyBytes, ImageFormat::RGBA8(), 512, 1, 1, 1, ImageFormat::RGBA16F());
    fb->set(Framebuffer::COLOR0, emissiveTex);
    rd->push2D(fb); {
        rd->setColorClearValue(Color3::black());
        rd->clear();
    } rd->pop2D();
    emissiveTex->generateMipMaps();
    spec.setEmissive(emissiveTex);
    //spec.setBump(System::findDataFile("material/10538-bump.jpg"));
    return UniversalMaterial::create(spec);
}

App::App(const GApp::Settings& settings) : GApp(settings) {
    renderDevice->setColorClearValue(Color3::white());
}

void App::onCleanup() {
	m_rtAudio.stopStream();
	if (m_rtAudio.isStreamOpen()) {
		m_rtAudio.closeStream();
	}
}

int audioCallback( void * outputBuffer, void * inputBuffer, unsigned int numFrames,
            double streamTime, RtAudioStreamStatus status, void * data ) {
	size_t numBytes = numFrames * sizeof(Sample);
	Array<Sample> output;
	output.resize(numFrames);
	for (Sample& s : output) {
		s = 0.0;
	}
	Synthesizer::global->synthesize(output);
	memcpy(outputBuffer, output.getCArray(), numBytes);

	// TODO: remove g_currentAudioBuffer... our mutex scheme handles everything just fine
	// Copy input buffer, handle race condition by not caring about it.
	memcpy(g_currentAudioBuffer.getCArray(), inputBuffer, numBytes);

    // For debugging
    /*
    for (int i = 0; i < numFrames; ++i) {
        int sampleIndex = (numFrames*g_sampleWindowIndex + i);
        //debugPrintf("sampleIndex: %d\n", sampleIndex);
        g_currentAudioBuffer[i] = sin(2*pif()*sampleIndex*(1.0/511.0));
    }*/

	s_app->m_rawAudioMutex.lock(); {
		++g_sampleWindowIndex;
		s_app->m_cpuRawAudioData.appendPOD(g_currentAudioBuffer);
        //debugPrintf("g_currentAudioBuffer.size() = %d\nm_cpuRawAudioData.size() = %d\n", g_currentAudioBuffer.size(), s_app->m_cpuRawAudioData.size());
	}s_app->m_rawAudioMutex.unlock();
	return 0;
}



void App::initializeAudio() {

  unsigned int bufferByteCount = 0;
  unsigned int bufferFrameCount = 512;
    
  // Check for audio devices
  if( m_rtAudio.getDeviceCount() < 1 ) {
    // None :(
    debugPrintf("No audio devices found!\n");
    exit( 1 );
  }
  RtAudio::DeviceInfo info;

  unsigned int devices = m_rtAudio.getDeviceCount();


  // Let RtAudio print messages to stderr.
  m_rtAudio.showWarnings( true );

  // Set input and output parameters
  RtAudio::StreamParameters iParams, oParams;
  iParams.deviceId = m_rtAudio.getDefaultInputDevice();
  iParams.nChannels = m_audioSettings.numChannels;
  iParams.firstChannel = 0;
  oParams.deviceId = m_rtAudio.getDefaultOutputDevice();
  oParams.nChannels = m_audioSettings.numChannels;
  oParams.firstChannel = 0;
    
  // Create stream options
  RtAudio::StreamOptions options;

  g_currentAudioBuffer.resize(bufferFrameCount);
  try {
    // Open a stream
    m_rtAudio.openStream( &oParams, &iParams, m_audioSettings.rtAudioFormat, m_audioSettings.sampleRate, &bufferFrameCount, &audioCallback, (void *)&bufferByteCount, &options );
  } catch( RtAudioError& e ) {
    // Failed to open stream
    std::cout << e.getMessage() << std::endl;
    exit( 1 );
  }
  g_currentAudioBuffer.resize(bufferFrameCount);
  m_rtAudio.startStream();

}

void App::onInit() {
    GApp::onInit();
    s_app = this;

    m_sonicSculptureFilename = "default.Soundscape.Any";

    m_freezeEverything = false;

    g_sampleWindowIndex = 0;
    m_lastSampleWindowProcessed = 0;
    m_initialTime = System::time();

    m_appMode = AppMode::DEFAULT;
    // TODO: Print instructions
    m_maxSavedTimeSlices = 128;
    m_lastInterestingEventTime = System::time();
    initializeAudio();

    m_rawAudioTexture = Texture::createEmpty("Raw Audio Texture", g_currentAudioBuffer.size(), 1, ImageFormat::R32F());
    m_frequencyAudioTexture = Texture::createEmpty("Frequency Audio Texture", g_currentAudioBuffer.size()/2, 1, ImageFormat::RG32F());

    m_smoothedRootMeanSquare = 0.0f;


    setFrameDuration(1.0f/30.0f);
    makeGUI();
    loadScene("Sculpting");
}

void App::makeGUI() {
    // Initialize the developer HUD (using the existing scene)
    createDeveloperHUD();
    debugWindow->setVisible(false);
    developerWindow->videoRecordDialog->setEnabled(true);
    debugPane->pack();


    debugWindow->pack();
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
    developerWindow->cameraControlWindow->setVisible(false);
    developerWindow->sceneEditorWindow->setVisible(false);
    developerWindow->setVisible(false);
    showRenderingStats = false;
    developerWindow->cameraControlWindow->moveTo(Point2(developerWindow->cameraControlWindow->rect().x0(), 0));


}


bool App::onEvent(const GEvent& e) {
    if (GApp::onEvent(e)) {
        return true;
    }

    // If you need to track individual UI events, manage them here.
    // Return true if you want to prevent other parts of the system
    // from observing this specific event.
    //
    // For example,
    // if ((e.type == GEventType::GUI_ACTION) && (e.gui.control == m_button)) { ... return true;}
    // if ((e.type == GEventType::KEY_DOWN) && (e.key.keysym.sym == GKey::TAB)) { ... return true; }

    return false;
}

void App::setAudioShaderArgs(Args& args) {
    m_rawAudioTexture->setShaderArgs(args, "rawAudio_", Sampler::video());
    m_frequencyAudioTexture->setShaderArgs(args, "frequencyAudio_", Sampler::video());
    m_fastMovingAverage.gpuData->setShaderArgs(args, "fastEWMAfreq_", Sampler::video());
    m_slowMovingAverage.gpuData->setShaderArgs(args, "slowEWMAfreq_", Sampler::video());
    m_glacialMovingAverage.gpuData->setShaderArgs(args, "glacialEWMAfreq_", Sampler::video());
}

void App::updateAudioData() {
    int sampleCount = g_currentAudioBuffer.size();
    int freqCount = sampleCount / 2;
    int numStoredTimeSlices = 0;
    int numNewWindows = 0;
    m_rawAudioMutex.lock(); {
        numStoredTimeSlices = m_cpuRawAudioData.size() / sampleCount;
        if (numStoredTimeSlices >= m_maxSavedTimeSlices) {
            int cutoffCount = numStoredTimeSlices - (m_maxSavedTimeSlices - 1);
            int newTotalSampleCount = sampleCount*(m_maxSavedTimeSlices - 1);

            for (int i = 0; i < newTotalSampleCount; ++i) {
                m_cpuRawAudioData[i] = m_cpuRawAudioData[i + sampleCount*cutoffCount];
            }
            m_cpuRawAudioData.resize(newTotalSampleCount);
            numStoredTimeSlices = m_cpuRawAudioData.size() / sampleCount;
        }
        m_cpuRawAudioSnapshot.fastClear();
        m_cpuRawAudioSnapshot.appendPOD(m_cpuRawAudioData);
        numNewWindows = min(g_sampleWindowIndex - m_lastSampleWindowProcessed, m_cpuRawAudioSnapshot.size() / sampleCount);
        m_lastSampleWindowProcessed = g_sampleWindowIndex;
    } m_rawAudioMutex.unlock();


	for (int n = numNewWindows - 1; n >= 0; --n) {
		float sumSquare = 0.0f;
		for (int i = m_cpuRawAudioSnapshot.size() - (sampleCount*(n+1)); i < m_cpuRawAudioSnapshot.size() - (sampleCount*n); ++i) {
			sumSquare += square(m_cpuRawAudioSnapshot[i]);
		}
		float rms = sqrt(sumSquare / sampleCount);
		m_smoothedRootMeanSquare = lerp(rms, m_smoothedRootMeanSquare, 0.8f);

		updateSonicSculpture(m_cpuRawAudioSnapshot.size() - (sampleCount*(n + 1)), sampleCount);
	}
		
	shared_ptr<CPUPixelTransferBuffer> ptb = CPUPixelTransferBuffer::fromData(sampleCount, numStoredTimeSlices, ImageFormat::R32F(), m_cpuRawAudioSnapshot.getCArray());
	m_rawAudioTexture->resize(sampleCount, numStoredTimeSlices);
	m_rawAudioTexture->update(ptb);

		
	
}

void App::handlePlayPulses() {
    for (int i = m_currentPlayPulses.size() - 1; i >= 0; --i) {
        int currentSampleIndex = (g_sampleWindowIndex * g_currentAudioBuffer.size());
        shared_ptr<SonicSculpturePiece> piece = m_currentPlayPulses[i].piece;
        int endIndex = m_currentPlayPulses[i].initialSample + (piece->size() * g_currentAudioBuffer.size());

        RenderDevice* rd = RenderDevice::current;
        static shared_ptr<Framebuffer> playPulseFB = Framebuffer::create("Play Pulse FB");
        shared_ptr<UniversalMaterial> material = piece->material();
        if (currentSampleIndex >= endIndex) {
            
            playPulseFB->set(Framebuffer::COLOR0, material->emissive().texture());
            rd->push2D(playPulseFB); {
                rd->setColorClearValue(Color3::black());
                rd->clear();
            } rd->pop2D();
            material->emissive().texture()->generateMipMaps();
            m_currentPlayPulses.remove(i);
            continue;
        }
        float alpha = float(currentSampleIndex - m_currentPlayPulses[i].initialSample) / (endIndex - m_currentPlayPulses[i].initialSample);
        
       
        playPulseFB->set(Framebuffer::COLOR0, material->emissive().texture());
        rd->push2D(playPulseFB); {
            Args args;
            args.setUniform("pulsePos", alpha * playPulseFB->width());
            args.setRect(rd->viewport());
            LAUNCH_SHADER("playPulse.pix", args);
        } rd->pop2D();
        material->emissive().texture()->generateMipMaps();
    }
}

void App::generatePlayPulse(shared_ptr<SonicSculpturePiece> piece) {
    PlayPulse pp;
    pp.initialSample = g_currentAudioBuffer.size() * g_sampleWindowIndex;
    pp.piece = piece;
    m_currentPlayPulses.append(pp);
    Synthesizer::global->queueSound(piece->getBaseAudioSample());
}

void App::saveSoundscape() const {
    Any a(Any::TABLE, "SoundScape");

    Array<Any> sonicSculptureArray;
    for (int i = 0; i < m_sonicSculpturePieces.size(); ++i) {
        sonicSculptureArray.append(m_sonicSculpturePieces[i]->toAny(format("sonicSculpture%d.bin", i)));
    }

    a["sonicSculptures"] = sonicSculptureArray;
    
    a.save(m_sonicSculptureFilename);
}

void App::loadSoundscape() {
    Any a;
    a.load(System::findDataFile(m_sonicSculptureFilename));
    AnyTableReader r("SoundScape", a);
    Array<Any> sonicSculptureArray;
    r.getIfPresent("sonicSculptures", sonicSculptureArray);
    r.verifyDone();
    m_sonicSculpturePieces.fastClear();
    for (int i = 0; i < sonicSculptureArray.size(); ++i) {
        shared_ptr<UniversalMaterial> material = getSonicSculptureMaterial(m_sonicSculpturePieces.size());
        m_sonicSculpturePieces.append(SonicSculpturePiece::create(material, sonicSculptureArray[i]));
    }
    
}


void App::onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& allSurfaces) {
    // This implementation is equivalent to the default GApp's. It is repeated here to make it
    // easy to modify rendering. If you don't require custom rendering, just delete this
    // method from your application and rely on the base class.

    if (!scene()) {
        if ((submitToDisplayMode() == SubmitToDisplayMode::MAXIMIZE_THROUGHPUT) && (!rd->swapBuffersAutomatically())) {
            swapBuffers();
        }
        rd->clear();
        rd->pushState(); {
            rd->setProjectionAndCameraMatrix(activeCamera()->projection(), activeCamera()->frame());
            drawDebugShapes();
        } rd->popState();
        return;
    }

    updateAudioData();

    if (System::time() - m_lastInterestingEventTime > 10.0) {
        if (Random::common().uniform() < 0.001f) {
            if (m_sonicSculpturePieces.size() > 0) {
                int index = Random::common().integer(0, m_sonicSculpturePieces.size() - 1);
                generatePlayPulse(m_sonicSculpturePieces[index]);
                m_lastInterestingEventTime = System::time();
            }
        }
    }
    handlePlayPulses();

    GBuffer::Specification gbufferSpec = m_gbufferSpecification;
    extendGBufferSpecification(gbufferSpec);
    m_gbuffer->setSpecification(gbufferSpec);
    m_gbuffer->resize(m_framebuffer->width(), m_framebuffer->height());
    m_gbuffer->prepare(rd, activeCamera(), 0, -(float)previousSimTimeStep(), m_settings.depthGuardBandThickness, m_settings.colorGuardBandThickness);

    m_renderer->render(rd, m_framebuffer, m_depthPeelFramebuffer, scene()->lightingEnvironment(), m_gbuffer, allSurfaces);

    // Debug visualizations and post-process effects
    rd->pushState(m_framebuffer); {
        // Call to make the App show the output of debugDraw(...)
        rd->setProjectionAndCameraMatrix(activeCamera()->projection(), activeCamera()->frame());

        for (auto& piece : m_sonicSculpturePieces) {
            piece->render(rd, scene()->lightingEnvironment());
        }
        if (notNull(m_currentSonicSculpturePiece)) {
            m_currentSonicSculpturePiece->render(rd, scene()->lightingEnvironment());
        }
        for (int i = m_playPlanes.size() - 1; i >= 0; --i) {
	        const PlayPlane& pp = m_playPlanes[i];
	        if (pp.endWindowIndex < g_sampleWindowIndex) {
	            m_playPlanes.remove(i);
	        }

	        Point3 point = pp.origin + (pp.direction * METERS_PER_SAMPLE_WINDOW * (g_sampleWindowIndex-pp.beginWindowIndex));

	        Color4 solidColor(1.0f, .02f, .03f, .15f);
	        //	  Plane plane(point, pp.direction);
	        //	  Draw::plane(plane, rd, solidColor, Color4::clear());

	        CFrame planeFrame = pp.frame;
	        planeFrame.translation = point;
	  
	        Box b(Vector3(-100.0f, -100.0f, 0.0f), Vector3(100.0f, 100.0f, 0.1f), planeFrame);
	        Draw::box(b, rd, solidColor, Color4::clear());
        }


        drawDebugShapes();
        const shared_ptr<Entity>& selectedEntity = (notNull(developerWindow) && notNull(developerWindow->sceneEditorWindow)) ? developerWindow->sceneEditorWindow->selectedEntity() : shared_ptr<Entity>();
        scene()->visualize(rd, selectedEntity, allSurfaces, sceneVisualizationSettings(), activeCamera());

        // Post-process special effects
        m_depthOfField->apply(rd, m_framebuffer->texture(0), m_framebuffer->texture(Framebuffer::DEPTH), activeCamera(), m_settings.depthGuardBandThickness - m_settings.colorGuardBandThickness);

        m_motionBlur->apply(rd, m_framebuffer->texture(0), m_gbuffer->texture(GBuffer::Field::SS_EXPRESSIVE_MOTION),
            m_framebuffer->texture(Framebuffer::DEPTH), activeCamera(),
            m_settings.depthGuardBandThickness - m_settings.colorGuardBandThickness);
    } rd->popState();

    rd->push2D(m_framebuffer); {
        rd->setBlendFunc(RenderDevice::BLEND_SRC_ALPHA, RenderDevice::BLEND_ONE_MINUS_SRC_ALPHA);
        Array<SoundInstance> soundInstances;
        Synthesizer::global->getSoundInstances(soundInstances);
        int xOffset = 10;
        Vector2 dim(256,100);
        for (int i = 0; i < soundInstances.size(); ++i) {
            int yOffset = rd->height() - 120 - (120 * i);
            Draw::rect2D(Rect2D::xywh(Point2(xOffset, yOffset), dim), rd, Color3::white(), soundInstances[i].displayTexture());
            float playheadAlpha = ((float)soundInstances[i].currentPosition) / soundInstances[i].audioSample->buffer.size();
            float playheadX = xOffset + (playheadAlpha * dim.x);
            Draw::rect2D(Rect2D::xywh(Point2(playheadX, yOffset), Vector2(1, dim.y)), rd, Color3::yellow());
        }


	
	static shared_ptr<GFont> font = GFont::fromFile(System::findDataFile("arial.fnt"));
	float time = System::time() - m_initialTime;
	if (time < 10.0f) {
	  float fontAlpha = time < 9.0f ? 1.0f : 10.0f - time;
	  font->draw2D(rd, "Press Space to Sculpt", Vector2(rd->width()/2, rd->height()-100.0f), 30.0f, Color4(Color3::black(), fontAlpha), Color4(Color3::white()*0.6f, fontAlpha), GFont::XALIGN_CENTER);
	}
    } rd->pop2D();

    if ((submitToDisplayMode() == SubmitToDisplayMode::MAXIMIZE_THROUGHPUT) && (!renderDevice->swapBuffersAutomatically())) {
        // We're about to render to the actual back buffer, so swap the buffers now.
        // This call also allows the screenshot and video recording to capture the
        // previous frame just before it is displayed.
        swapBuffers();
    }

    // Clear the entire screen (needed even though we'll render over it, since
    // AFR uses clear() to detect that the buffer is not re-used.)
    rd->clear();

    // Perform gamma correction, bloom, and SSAA, and write to the native window frame buffer
    m_film->exposeAndRender(rd, activeCamera()->filmSettings(), m_framebuffer->texture(0));
}


void App::playSculpture(const Ray& playRay) {
  int maxDistance = 0;
  int startIndex = g_sampleWindowIndex;
  for (auto piece : m_sonicSculpturePieces) {
    const shared_ptr<AudioSample>& sample = piece->getAudioSampleFromRay(playRay);
    if (sample->buffer.size() > 0) {
        maxDistance = max(maxDistance, (int)ceil(sample->buffer.size() / 512.0f));
        Synthesizer::global->queueSound(sample);
        m_lastInterestingEventTime = System::time();
    }
  }
  if (maxDistance > 0) {
    PlayPlane pp;
    pp.direction = playRay.direction();
    pp.origin = playRay.origin();
    pp.beginWindowIndex = startIndex;
    pp.endWindowIndex = startIndex + maxDistance;

    Vector3 zAxis = -playRay.direction();
	Vector3 xAxis;
	if (abs(zAxis.dot(Vector3::unitY())) < 0.9f) {
	  xAxis = zAxis.unitCross(Vector3::unitY());
	} else {
	  xAxis = zAxis.unitCross(Vector3::unitX());
	}
	Vector3 yAxis = zAxis.cross(xAxis);

	pp.frame = CFrame(Matrix3::fromColumns(xAxis, yAxis, zAxis), pp.origin);


    m_playPlanes.append(pp);
  }
}

void App::onUserInput(UserInput* ui) {
    GApp::onUserInput(ui);
    if (ui->keyDown(GKey::SPACE)) {
        m_appMode = AppMode::MAKING_SCULPTURE;
    } else {
        m_appMode = AppMode::DEFAULT;
    }

    if (ui->keyDown(GKey::LCTRL) || ui->keyDown(GKey::RCTRL)) {
        if (ui->keyPressed(GKey('s'))) {
            saveSoundscape();
        }
        if (ui->keyPressed(GKey('l'))) {
            loadSoundscape();
        }
    }

	if (ui->keyPressed(GKey::RETURN)) {


        Point3 origin = Point3(0.0, 0.0, 0.0);
        Vector3 direction = activeCamera()->frame().rightVector();
        Ray playRay(origin, direction);
        float minValue = finf();
        for (auto& sculpture : m_sonicSculpturePieces) {
            minValue = min(minValue, sculpture->minValueAlongRay(playRay));
        }

        playRay.set(origin + direction*minValue, direction);
		playSculpture(playRay);
	}

    if (ui->keyPressed(GKey('f'))) {
        m_freezeEverything = !m_freezeEverything;
        for (shared_ptr<SonicSculpturePiece>& piece : m_sonicSculpturePieces) {
            piece->setFrozen(m_freezeEverything);
        }
    }

    // Hack for playing pieces
    for (int i = 0; i < 8; ++i) {
        if (m_sonicSculpturePieces.size() > i) {
            if (ui->keyPressed(GKey('1' + char(i)))) {
                generatePlayPulse(m_sonicSculpturePieces[i]);
                m_lastInterestingEventTime = System::time();
            }
        }
    }


    // Add key handling here based on the keys currently held or
    // ones that changed in the last frame.
}

void App::onGraphics2D(RenderDevice* rd, Array<Surface2D::Ref>& posed2D) {
    // Render 2D objects like Widgets.  These do not receive tone mapping or gamma correction
    Surface2D::sortAndRender(rd, posed2D);
}


void App::updateSonicSculpture(int audioSampleOffset, int audioSampleCount) {
	float delta = 0.1f;
	float radius = sqrt(m_smoothedRootMeanSquare) * 1.0f;
	CFrame frame = activeCamera()->frame();
	// Offset a bit forward
	frame.translation += activeCamera()->frame().lookVector() * 0.2f;
	// Offset 3 inches down to the mouth opening
	frame.translation += activeCamera()->frame().upVector() * -0.0762f;
	if (m_appMode == AppMode::DEFAULT) {
		if (notNull(m_currentSonicSculpturePiece)) {
			if (m_currentSonicSculpturePiece->size() > 0) {
                m_sonicSculpturePieces.append(m_currentSonicSculpturePiece);
                m_lastInterestingEventTime = System::time();
			}
			m_currentSonicSculpturePiece = shared_ptr<SonicSculpturePiece>();
		}
	} else if (m_appMode == AppMode::MAKING_SCULPTURE) {
		if (isNull(m_currentSonicSculpturePiece)) {
            shared_ptr<UniversalMaterial> material = getSonicSculptureMaterial(m_sonicSculpturePieces.size());
			m_currentSonicSculpturePiece = SonicSculpturePiece::create(material);
		}
        m_lastInterestingEventTime = System::time();
		// TODO: eliminate some of the redundant copies
		Array<float> samples;
		samples.resize(audioSampleCount);
		for (int i = 0; i < audioSampleCount; ++i) {
			samples[i] = m_cpuRawAudioSnapshot[i + audioSampleOffset];
		}
		m_currentSonicSculpturePiece->insert(frame, radius, delta, samples);
	}
}

void App::onSimulation(RealTime rdt, SimTime sdt, SimTime idt) {
    GApp::onSimulation(rdt, sdt, idt);

    for (int i = 0; i < m_sonicSculpturePieces.size(); ++i) {
        m_sonicSculpturePieces[i]->onSimulation(rdt, sdt, idt);
    }
    if (notNull(m_currentSonicSculpturePiece)) {
        m_currentSonicSculpturePiece->onSimulation(rdt, sdt, idt);
    }
    // See 
    Synthesizer::global->flushHackQueue();

    // Example GUI dynamic layout code.  Resize the debugWindow to fill
    // the screen horizontally.
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}
