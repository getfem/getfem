#ifndef STREAM_REDIRECT_H
#define STREAM_REDIRECT_H

#include <iostream>
#include <streambuf>
#include <string>

//////////////////////////
// For cout redirection //
//////////////////////////

class ScilabStream : public std::basic_streambuf<char>
{
public:
  ScilabStream(std::ostream &stream) : m_stream(stream)
  {
    m_old_buf = stream.rdbuf();
    stream.rdbuf(this);
  }
  ~ScilabStream()
  {
    // output anything that is left
    if (!m_string.empty())
      sciprint("getfem: %s\n",m_string.c_str());

    m_stream.rdbuf(m_old_buf);
  }

protected:
  virtual int_type overflow(int_type v)
  {
    if (v == '\n')
      {
	sciprint("getfem: %s\n",m_string.c_str());
	m_string.clear();
      }
    else
      m_string.push_back(v);
    
    return v;
  }
  
  virtual std::streamsize xsputn(const char *p, std::streamsize n) 
  {
    m_string.append(p, p + n);
    
    int pos = 0;
    while (pos != std::string::npos)
      {
	pos = m_string.find('\n');
	if (pos != std::string::npos)
	  {
	    std::string tmp(m_string.begin(), m_string.begin() + pos);
	    sciprint("getfem: %s\n",tmp.c_str());
	    m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
	  }
      }
    
    return n;
  }
  
private:
  std::ostream   &m_stream;
  std::streambuf *m_old_buf;
  std::string     m_string;
};
#endif
