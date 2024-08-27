read_biolector <- function(
  file,
  measureType = "AMP",
  blVersion = 1
)
{
  if(blVersion==3) {
    dat = read_biolector3(file, measureType=measureType)
  }
  else {
    dat = read_biolector1(file, measureType=measureType)
  }
  return(dat)
}