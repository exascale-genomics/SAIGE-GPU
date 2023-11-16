# Shrink Wrap
A std::streambuf wrapper for compression formats.

## XZ streambuf with std::istream
```c++
std::array<char, 1024> buf;
shrinkwrap::xz::ibuf sbuf("file.xz");
std::istream is(&sbuf);
is.seekg(-1024, std::ios::end);
while (is)
{
  is.read(buf.data(), buf.size());
  std::cout.write(buf.data(), is.gcount());
}
```
## XZ streambuf with std::istreambuf_iterator
```c++
shrinkwrap::xz::ibuf sbuf("file.xz");
for (std::istreambuf_iterator<char> it(&sbuf); it != std::istreambuf_iterator<char>{}; ++it)
  std::cout.put(*it);
```

## XZ input stream 
```c++
std::array<char, 1024> buf;
shrinkwrap::xz::ostream is("file.xz");
while (is)
{
  is.read(buf.data(), buf.size());
  std::cout.write(buf.data(), is.gcount());
}
```

## XZ output stream 
```c++
std::vector<char> buf(1024 * 1024);
shrinkwrap::xz::ostream os("file.xz");
while (std::cin)
{
  std::cin.read(buf.data(), buf.size());
  os.write(buf.data(), buf.gcount());
  os.flush(); // flush() creates block boundary.
}
```

## BGZF (Blocked GNU Zip Format)  
```c++
std::array<char, 1024> buf;
shrinkwrap::bgzf::istream is("file.xz");
is.read(buf.data(), buf.size());

// (gzip_block_position << 16) | relative_uncompressed_offset
auto virtual_offset = is.tellg();
is.seekg(virtual_offset);
```

## Generic input stream
Generic istream detects file format.
```c++
std::array<char, 1024> buf;
shrinkwrap::istream is("file");
while (is)
{
  is.read(buf.data(), buf.size());
  std::cout.write(buf.data(), is.gcount());
}
```

## Caveats
* Does not support files with concatenated xz streams.
