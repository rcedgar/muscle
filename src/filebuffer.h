/////////////////////////////////////////////////////////////////
// filebuffer.h
//
// Buffered file reading.
/////////////////////////////////////////////////////////////////


#ifndef FILEBUFFER_H
#define FILEBUFFER_H

#include <string>
#include <fstream>
#include <iostream>
#include <assert.h>

using namespace std;

const int BufferSize = 1000;

/////////////////////////////////////////////////////////////////
// FileBuffer
//
// Class for buffering file reading.
/////////////////////////////////////////////////////////////////

class FileBuffer {
  ifstream file;
  char buffer[BufferSize];
  int currPos;
  int size;
  bool isEOF;
  bool isValid;
  bool canUnget;

 public:

  // Some common routines

  FileBuffer (const char *filename) : file (filename), currPos (0), size (0), isEOF (false), isValid (!file.fail()), canUnget (false){}
  ~FileBuffer (){ close(); }
  bool fail () const { return !isValid; }
  bool eof () const { return (!isValid || isEOF); }
  void close(){ file.close(); isValid = false; }

  /////////////////////////////////////////////////////////////////
  // FileBuffer::Get()
  //
  // Retrieve a character from the file buffer.  Returns true if
  // and only if a character is read.
  /////////////////////////////////////////////////////////////////

  bool Get (char &ch){

    // check to make sure that there's more stuff in the file
    if (!isValid || isEOF) return false;

    // if the buffer is empty, it's time to reload it
    if (currPos == size){
      file.read (buffer, BufferSize);
      size = int(file.gcount());
      isEOF = (size == 0);
      currPos = 0;
      if (isEOF) return false;
    }

    // store the read character
    ch = buffer[currPos++];
    canUnget = true;
    return true;
  }

  /////////////////////////////////////////////////////////////////
  // FileBuffer::UnGet()
  //
  // Unretrieve the most recently read character from the file
  // buffer.  Note that this allows only a one-level undo.
  /////////////////////////////////////////////////////////////////

  void UnGet (){
    assert (canUnget);
    assert (isValid);
    assert (currPos > 0);
    currPos--;
    assert (currPos < size);
    isEOF = false;
    canUnget = false;
  }

  /////////////////////////////////////////////////////////////////
  // FileBuffer::GetLine()
  //
  // Retrieve characters of text until a newline character is
  // encountered.  Terminates properly on end-of-file condition.
  /////////////////////////////////////////////////////////////////

  void GetLine (string &s){
    char ch;
    s = "";
    while (Get (ch) && ch != '\n')
      s += ch;
  }

};

#endif
