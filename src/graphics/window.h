#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include "shaders.h"

#include "../sources/source.h"
#include "../detectors/detector.h"
#include "../fields/fields.h"

static GLFWwindow* window;
static GLuint VertexArrayID, TextureID, TexUniformID, programID, MatrixID;
static GLuint vertexbuffers[2];

static glm::mat4 Projection, View, Model, mvp;

int openWindow(const char * vertex_file_path, const char * fragment_file_path);
int openWindowSized(const char * vertex_file_path, const char * fragment_file_path, int size[2]);

int setupMatrix(ParticleSource* source, FieldStructure* field, ParticleDetector* detector);
int setupMatrix(double max_beta, FieldStructure* field);
int setupBuffers(double pos[], bool running[], long N);
int updateBuffers(double pos[], bool running[], long N);
int draw(int N, int name_length, char* name, bool wait = false);

int setupBuffersTexRectangle(float* tbo_data, int pixels[2]);
int updateBuffersTexRectangle(float* tbo_data, int pixels[2]);
int drawTexRectangle();
