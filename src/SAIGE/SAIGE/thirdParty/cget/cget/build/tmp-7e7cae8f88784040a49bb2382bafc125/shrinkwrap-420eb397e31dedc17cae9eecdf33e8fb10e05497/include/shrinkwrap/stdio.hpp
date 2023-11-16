#ifndef SHRINKWRAP_STDIO_HPP
#define SHRINKWRAP_STDIO_HPP

#include <streambuf>
#include <cstdio>
#include <vector>
#include <iostream>
#include <cassert>

namespace shrinkwrap
{
  namespace stdio
  {
    class filebuf : public std::streambuf
    {
    private:
      std::FILE* fp_;
    public:
      typedef filebuf self_type;

      filebuf(std::FILE* fp) :
        fp_(fp)
      { }

      filebuf(self_type&& src) :
        std::streambuf(std::move(src)),
        fp_(src.fp_)
      {
        src.fp_ = nullptr;
      }

      self_type& operator=(self_type&& src)
      {
        if (this != &src)
        {
          std::streambuf::operator=(std::move(src));
          this->destroy();
          std::swap(fp_, src.fp_);
        }
        return *this;
      }

      void swap(self_type& src)
      {
        std::streambuf::swap(src);
        std::swap(fp_, src.fp_);
      }

      ~filebuf()
      {
        this->destroy();
      }

      void destroy()
      {
        if (fp_)
        {
          std::fclose(fp_);
          fp_ = nullptr;
        }
      }

      bool is_open() const
      {
        return fp_ != nullptr;
      }
    protected:
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // buffer management and positioning virtual methods
      virtual std::streambuf::pos_type seekoff(std::streambuf::off_type off, std::ios_base::seekdir way, std::ios_base::openmode which)
      {
        int origin;
        if (way == std::ios::cur)
          origin = SEEK_CUR;
        else if (way == std::ios::end)
          origin = SEEK_END;
        else
          origin = SEEK_SET;

        if (fseek(fp_, off, origin))
          return pos_type(off_type(-1));
        return pos_type(off_type(std::ftell(fp_)));
      }

      virtual std::streambuf::pos_type seekpos(std::streambuf::pos_type pos, std::ios_base::openmode which)
      {
        if (fseek(fp_, pos, SEEK_SET))
          return pos_type(off_type(-1));
        assert(pos == pos_type(off_type(std::ftell(fp_))));
        return pos;
      }

      virtual int sync()
      {
        return std::fflush(fp_);
      }
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // input virtual methods
      virtual std::streamsize xsgetn(char* data, std::streamsize sz)
      {
        std::streamsize ret = std::fread(data, 1, sz, fp_);

        if (std::ferror(fp_))
        {
          // TODO: throw exception
        }

        return ret;
      }

      virtual int underflow() // This supports istream::peek()
      {
        int ret = std::getc(fp_);
        return std::ungetc(ret, fp_);
      }

      virtual int uflow() // This supports istream::get()
      {
        return std::getc(fp_);
      }

      virtual int pbackfail(int c)
      {
        return std::ungetc(c, fp_); // Does not support istream::unget(), which is indicated by c equaling EOF
      }
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // output virtual methods
      virtual std::streamsize xsputn(const char* data, std::streamsize sz)
      {
        std::streamsize ret = std::fwrite(data, 1, sz, fp_);

        if (ret != sz)
        {
          // TODO throw exception
        }

        return ret;
      }

      virtual int overflow(int c)
      {
        if (!traits_type::eq_int_type(c, traits_type::eof()))
          return std::putc(c, fp_);
        return traits_type::not_eof(c);
      }
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    };

    class istream : public std::istream
    {
    public:
      istream(const std::string& file_path) :
        std::istream(&sbuf_),
        sbuf_(std::fopen(file_path.c_str(), "rb"))
      {
        if (!sbuf_.is_open())
          setstate(rdstate() | std::ios::badbit);
      }

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
      istream(istream&& src)
        :
        std::istream(&sbuf_),
        sbuf_(std::move(src.sbuf_))
      {
      }

      istream& operator=(istream&& src)
      {
        if (&src != this)
        {
          std::istream::operator=(std::move(src));
          sbuf_ = std::move(src.sbuf_);
        }
        return *this;
      }
#endif
    private:
      ::shrinkwrap::stdio::filebuf sbuf_;
    };



    class ostream : public std::ostream
    {
    public:
      ostream(const std::string& file_path)
        :
        std::ostream(&sbuf_),
        sbuf_(std::fopen(file_path.c_str(), "wb"))
      {
        if (!sbuf_.is_open())
          setstate(rdstate() | std::ios::badbit);
      }

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
      ostream(ostream&& src)
        :
        std::ostream(&sbuf_),
        sbuf_(std::move(src.sbuf_))
      {
      }

      ostream& operator=(ostream&& src)
      {
        if (&src != this)
        {
          std::ostream::operator=(std::move(src));
          sbuf_ = std::move(src.sbuf_);
        }
        return *this;
      }
#endif
    private:
      ::shrinkwrap::stdio::filebuf sbuf_;
    };
  }
}

#endif //SHRINKWRAP_STDIO_HPP