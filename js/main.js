/*
<!--[if lt IE 9]>
<sricpt src="http://html5shiv.googlecode.com/svn/trunk/html5.js"></script>
<![endif]-->
*/



// block IE (versions before 9.0)

{
    var userAgent = navigator.userAgent.toLowerCase();
    var browser = {
        version: (userAgent.match(/(?:firefox|opera|safari|chrome|msie)[\/: ]([\d.]+)/))[1],
        safari: /version.+safari/.test(userAgent),
        chrome: /chrome/.test(userAgent),
        firefox: /firefox/.test(userAgent),
        ie: /msie/.test(userAgent),
        opera: /opera/.test(userAgent)
    }
    if (browser.ie && browser.version < 9) {
        alert("Sorry, this site appears to be too cool for your web browser(>_<)\nPlease use one of these recommended browsers instead:\n\nGoogle Chrome\nApple Safari\nMozilla Firefox\nMicrosoft Internet Explorer 9.0 or later");
        window.location="https://www.google.com/chrome";
    }
}


//Add this to handle PDF links more gracefully:
document.querySelectorAll('a[href$=".pdf"]').forEach(link => {
  link.setAttribute('target', '_blank');
  link.setAttribute('rel', 'noopener noreferrer');
  link.insertAdjacentHTML('beforeend', ' <small>(PDF)</small>');
});