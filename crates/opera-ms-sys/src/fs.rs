//! Safe wrapper around `statvfs`.

use std::ffi::CString;
use std::io;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

/// The bytes available to an unprivileged user on the file system of
/// `path`.
pub fn free_space(path: &Path) -> io::Result<u64> {
    let c_path = CString::new(path.as_os_str().as_bytes())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "path contains a NUL byte"))?;
    let mut st = std::mem::MaybeUninit::<libc::statvfs>::uninit();
    // SAFETY: `c_path` is a valid NUL-terminated string and `st` points to
    // writable memory of the right size; `statvfs` only writes it.
    let rc = unsafe { libc::statvfs(c_path.as_ptr(), st.as_mut_ptr()) };
    if rc != 0 {
        return Err(io::Error::last_os_error());
    }
    // SAFETY: `statvfs` returned 0, so it filled `st`.
    let st = unsafe { st.assume_init() };
    // the field types differ between platforms (u64 on Linux, u32 on macOS)
    #[allow(clippy::useless_conversion)]
    Ok(u64::from(st.f_bavail).saturating_mul(u64::from(st.f_frsize)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn free_space_of_a_directory() {
        let dir = std::env::temp_dir();
        assert!(free_space(&dir).unwrap() > 0);
        assert!(free_space(&dir.join("opera-ms-no-such-dir")).is_err());
    }
}
